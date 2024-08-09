%% path continuation - dpa.checkResidual
%
%   Institute of Dynamics and Vibration Research
%   Leibniz University Hannover
%   04.04.2022 - Alwin Förster
%
function [dpaPoints] = checkResidual(fun,dpaPoints,oih)
    %% check dpaResidual
    n = oih.path.nAll;
    nv = numel(oih.path.varAll(:,1));
    isZero = 0;
    if n==2
        rnm1 = oih.opt.dpaResidual(oih.path.varAll(:,1),oih.path.lAll(1),oih.opt.g0);
        rn = oih.opt.dpaResidual(oih.path.varAll(:,2),oih.path.lAll(2),oih.opt.g0);
        if sign(rnm1*rn)<0
            isZero = 2;
        end
    elseif n>2
        rnm2 = oih.opt.dpaResidual(oih.path.varAll(:,n-2),oih.path.lAll(n-2),oih.opt.g0);
        rnm1 = oih.opt.dpaResidual(oih.path.varAll(:,n-1),oih.path.lAll(n-1),oih.opt.g0);
        rn = oih.opt.dpaResidual(oih.path.varAll(:,n),oih.path.lAll(n),oih.opt.g0);
        if sign(rn*rnm1)<0
            isZero = 2;
        elseif abs(rn)>abs(rnm1) && abs(rnm1)<abs(rnm2) && sign(rn*rnm1)>0 && sign(rnm1*rnm2)>0
            isZero = 1;
        end
    end
    %% find solution
    if isZero
        %% initial guess
        warning off;
        switch isZero
            case 1
                xl = [oih.path.varAll(:,n+[-2,-1,0]);oih.path.lAll(n+[-2,-1,0])];
                sl = oih.path.sAll(n+[-2,-1,0]);
                pl = polyfit(sl,[rnm2,rnm1,rn],2);
                dpl = polyder(pl,1);
                s0 = -dpl(2)/dpl(1);
                px = poly.fitn(sl,xl,2,true);
                x0 = poly.valn(px,s0,nv+1);
            case 2
                xl = [oih.path.varAll(:,n+[-1,0]);oih.path.lAll(n+[-1,0])];
                pl = poly.fitn([rnm1,rn],xl,1);
                x0 = poly.valn(pl,0,nv+1);
            otherwise
                error('Unknown error.');
        end
        warning on;
        %% solve
        isSolution = false;
        index = inf;
        Rres = @(x) dpa.mergeResiduals(oih,fun,oih.opt.dpaResidual,x,oih.opt.g0);
        [xr,~,exitflag] = oih.solver.numJac(Rres,x0);
        if exitflag>0
            switch isZero
                case 1
                    xnm2 = [oih.path.varAll(:,n-2);oih.path.lAll(n-2)];
                    xnm1 = [oih.path.varAll(:,n-1);oih.path.lAll(n-1)];
                    xn = [oih.path.varAll(:,n);oih.path.lAll(n)];
                    normD1 = norm(xn-xnm1);
                    normD2 = norm(xnm1-xnm2);
                    normNm2 = norm(xr-xnm2);
                    normNm1 = norm(xr-xnm1);
                    normN = norm(xr-xn);
                    if normD1>normNm1 && normD1>normN
                        isSolution = true;
                        index = n;
                    elseif normD2>normNm2 && normD2>normNm1
                        isSolution = true;
                        index = n-1;
                    end
                case 2
                    xnm1 = [oih.path.varAll(:,n-1);oih.path.lAll(n-1)];
                    xn = [oih.path.varAll(:,n);oih.path.lAll(n)];
                    normD = norm(xn-xnm1);
                    normNm1 = norm(xr-xnm1);
                    normN = norm(xr-xn);
                    if normD>normNm1 && normD>normN
                        isSolution = true;
                        index = n;
                    end
                otherwise
                    error('Unknown error.');
            end
        end
        if isSolution
            varAll = oih.path.varAll;
            lAll = oih.path.lAll;
            varAll = [varAll(:,1:(index-1)),xr(1:nv),varAll(:,index:end)];
            lAll = [lAll(1:(index-1)),xr(nv+1),lAll(index:end)];
            sAll = cumsum([0,sqrt(sum(([varAll(:,2:end);lAll(2:end)]-[varAll(:,1:(end-1));lAll(1:(end-1))]).^2,1))]);
            oih.path.varAll = varAll;
            oih.path.lAll = lAll;
            oih.path.sAll = sAll;
            dpaPoints = [dpaPoints,index];
        end
    end
end