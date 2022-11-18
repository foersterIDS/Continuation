%% path continuation - dpa.check_residual
%
%   Institute of Dynamics and Vibration Research
%   Leibniz University Hannover
%   04.04.2022 - Alwin Förster
%
function [Path,dpa_points] = check_residual(fun,dpa_points,Opt,Path,Solver)
    %% check dpa_residual
    n = numel(Path.l_all);
    nv = numel(Path.var_all(:,1));
    is_zero = 0;
    if n==2
        rnm1 = Opt.dpa_residual(Path.var_all(:,1),Path.l_all(1),Opt.g_0);
        rn = Opt.dpa_residual(Path.var_all(:,2),Path.l_all(2),Opt.g_0);
        if sign(rnm1*rn)<0
            is_zero = 2;
        end
    elseif n>2
        rnm2 = Opt.dpa_residual(Path.var_all(:,n-2),Path.l_all(n-2),Opt.g_0);
        rnm1 = Opt.dpa_residual(Path.var_all(:,n-1),Path.l_all(n-1),Opt.g_0);
        rn = Opt.dpa_residual(Path.var_all(:,n),Path.l_all(n),Opt.g_0);
        if sign(rn*rnm1)<0
            is_zero = 2;
        elseif abs(rn)>abs(rnm1) && abs(rnm1)<abs(rnm2) && sign(rn*rnm1)>0 && sign(rnm1*rnm2)>0
            is_zero = 1;
        end
    end
    %% find solution
    if is_zero
        %% initial guess
        warning off;
        switch is_zero
            case 1
                xl = [Path.var_all(:,n+[-2,-1,0]);Path.l_all(n+[-2,-1,0])];
                sl = Path.s_all(n+[-2,-1,0]);
                pl = polyfit(sl,[rnm2,rnm1,rn],2);
                dpl = polyder(pl,1);
                s0 = -dpl(2)/dpl(1);
                px = poly.fitn(sl,xl,2,true);
                x0 = poly.valn(px,s0,nv+1);
            case 2
                xl = [Path.var_all(:,n+[-1,0]);Path.l_all(n+[-1,0])];
                pl = poly.fitn([rnm1,rn],xl,1);
                x0 = poly.valn(pl,0,nv+1);
            otherwise
                error('Unknown error.');
        end
        warning on;
        %% solve
        is_solution = false;
        index = inf;
        Rres = @(x) dpa.merge_residuals(Opt,fun,Opt.dpa_residual,x,Opt.g_0);
        [xr,~,exitflag] = Solver.num_jac(Rres,x0);
        if exitflag>0
            switch is_zero
                case 1
                    xnm2 = [Path.var_all(:,n-2);Path.l_all(n-2)];
                    xnm1 = [Path.var_all(:,n-1);Path.l_all(n-1)];
                    xn = [Path.var_all(:,n);Path.l_all(n)];
                    norm_d1 = norm(xn-xnm1);
                    norm_d2 = norm(xnm1-xnm2);
                    norm_nm2 = norm(xr-xnm2);
                    norm_nm1 = norm(xr-xnm1);
                    norm_n = norm(xr-xn);
                    if norm_d1>norm_nm1 && norm_d1>norm_n
                        is_solution = true;
                        index = n;
                    elseif norm_d2>norm_nm2 && norm_d2>norm_nm1
                        is_solution = true;
                        index = n-1;
                    end
                case 2
                    xnm1 = [Path.var_all(:,n-1);Path.l_all(n-1)];
                    xn = [Path.var_all(:,n);Path.l_all(n)];
                    norm_d = norm(xn-xnm1);
                    norm_nm1 = norm(xr-xnm1);
                    norm_n = norm(xr-xn);
                    if norm_d>norm_nm1 && norm_d>norm_n
                        is_solution = true;
                        index = n;
                    end
                otherwise
                    error('Unknown error.');
            end
        end
        if is_solution
            var_all = Path.var_all;
            l_all = Path.l_all;
            var_all = [var_all(:,1:(index-1)),xr(1:nv),var_all(:,index:end)];
            l_all = [l_all(1:(index-1)),xr(nv+1),l_all(index:end)];
            s_all = cumsum([0,sqrt(sum(([var_all(:,2:end);l_all(2:end)]-[var_all(:,1:(end-1));l_all(1:(end-1))]).^2,1))]);
            Path.var_all = var_all;
            Path.l_all = l_all;
            Path.s_all = s_all;
            dpa_points = [dpa_points,index];
        end
    end
end