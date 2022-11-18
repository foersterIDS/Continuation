%% path continuation - dpa.checkResidual
%
%   Institute of Dynamics and Vibration Research
%   Leibniz University Hannover
%   04.04.2022 - Alwin Förster
%
function [Path,dpaPoints] = checkResidual(fun,dpaPoints,Opt,Path,Solver)
    %% check dpaResidual
    n = numel(Path.lAll);
    nv = numel(Path.varAll(:,1));
    isZero = 0;
    if n==2
        rnm1 = Opt.dpaResidual(Path.varAll(:,1),Path.lAll(1),Opt.g0);
        rn = Opt.dpaResidual(Path.varAll(:,2),Path.lAll(2),Opt.g0);
        if sign(rnm1*rn)<0
            isZero = 2;
        end
    elseif n>2
        rnm2 = Opt.dpaResidual(Path.varAll(:,n-2),Path.lAll(n-2),Opt.g0);
        rnm1 = Opt.dpaResidual(Path.varAll(:,n-1),Path.lAll(n-1),Opt.g0);
        rn = Opt.dpaResidual(Path.varAll(:,n),Path.lAll(n),Opt.g0);
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
                xl = [Path.varAll(:,n+[-2,-1,0]);Path.lAll(n+[-2,-1,0])];
                sl = Path.sAll(n+[-2,-1,0]);
                pl = polyfit(sl,[rnm2,rnm1,rn],2);
                dpl = polyder(pl,1);
                s0 = -dpl(2)/dpl(1);
                px = poly.fitn(sl,xl,2,true);
                x0 = poly.valn(px,s0,nv+1);
            case 2
                xl = [Path.varAll(:,n+[-1,0]);Path.lAll(n+[-1,0])];
                pl = poly.fitn([rnm1,rn],xl,1);
                x0 = poly.valn(pl,0,nv+1);
            otherwise
                error('Unknown error.');
        end
        warning on;
        %% solve
        isSolution = false;
        index = inf;
        Rres = @(x) dpa.mergeResiduals(Opt,fun,Opt.dpaResidual,x,Opt.g0);
        [xr,~,exitflag] = Solver.numJac(Rres,x0);
        if exitflag>0
            switch isZero
                case 1
                    xnm2 = [Path.varAll(:,n-2);Path.lAll(n-2)];
                    xnm1 = [Path.varAll(:,n-1);Path.lAll(n-1)];
                    xn = [Path.varAll(:,n);Path.lAll(n)];
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
                    xnm1 = [Path.varAll(:,n-1);Path.lAll(n-1)];
                    xn = [Path.varAll(:,n);Path.lAll(n)];
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
            varAll = Path.varAll;
            lAll = Path.lAll;
            varAll = [varAll(:,1:(index-1)),xr(1:nv),varAll(:,index:end)];
            lAll = [lAll(1:(index-1)),xr(nv+1),lAll(index:end)];
            sAll = cumsum([0,sqrt(sum(([varAll(:,2:end);lAll(2:end)]-[varAll(:,1:(end-1));lAll(1:(end-1))]).^2,1))]);
            Path.varAll = varAll;
            Path.lAll = lAll;
            Path.sAll = sAll;
            dpaPoints = [dpaPoints,index];
        end
    end
end