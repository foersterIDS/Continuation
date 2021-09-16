%% path continuation - predictor
%
%   Institute of Dynamics and Vibration Research
%   Leibniz University Hannover
%   08.05.2020 - Alwin Förster
%
function [vp,lp,fun_predictor,sp] = predictor(Path,ds,solver_jacobian,fun,res_corr,predictor_solver,Opt)
    %% get fun_predictor:
    if Opt.predictor.polynomial
        if length(Path.l_all)==1
            fun_predictor = @(s) predictor_initial(Path,s,Opt);
        else
            [nt,nf] = predictor_adaptive(Path,Opt);
            [fpt,Jpt] = predictor_taylor(Path,nt,nf);
            fun_predictor = @(s) fncHndToVal(s,fpt,Jpt);
        end
    elseif Opt.predictor.tangential
        if length(Path.l_all)==1
            fun_predictor = @(s) predictor_initial(Path,s,Opt);
        else
            fun_predictor = @(s) predictor_ode(Path,s,solver_jacobian,fun);
        end
    else
        error('predictor not set or of unknown type');
    end
    %% predictor_solver:
    if Opt.predictor_solver
        xi = [Path.var_all(:,end);Path.l_all(end)];
        fun_solve = @(s) merge_arle_pred(fun_predictor,res_corr,s,xi,ds);
        [sp,~,exitflag] = predictor_solver(fun_solve,ds);
        if exitflag<=0
            sp = ds;
        end
    else
        sp = ds;
    end
    %% get predictor:
    xip1 = fun_predictor(sp);
    %% correct_predictor:
    if Opt.correct_predictor && numel(Path.l_all)>1
        dxi = [Path.var_all(:,end);Path.l_all(end)]-[Path.var_all(:,end-1);Path.l_all(end-1)];
        dxip1 = xip1-[Path.var_all(:,end);Path.l_all(end)];
        if dot(dxip1,dxi)<0
            dxip1 = -dxip1;
            xip1 = [Path.var_all(:,end);Path.l_all(end)]+dxip1;
        end
    end
    %% make output:
    vp = xip1(1:end-1);
    lp = xip1(end);
end