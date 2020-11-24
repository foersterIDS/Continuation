%% path continuation - predictor
%
%   Institute of Dynamics and Vibration Research
%   Leibniz University Hannover
%   08.05.2020 - Alwin Förster
%
function [vp,lp,fun_predictor] = predictor(var_all,l_all,s_all,ds,solver_jacobian,fun,res_arle,predictor_solver,Opt)
    if Opt.predictor.polynomial
        if length(l_all)==1
            fun_predictor = @(s) predictor_initial(var_all,l_all,s,Opt);
        else
            [nt,nf] = predictor_adaptive(var_all,l_all,s_all,Opt);
            [fpt,Jpt] = predictor_taylor(var_all,l_all,s_all,nt,nf);
            fun_predictor = @(s) fncHndToVal(s,fpt,Jpt);
        end
    elseif Opt.predictor.tangential
        if length(l_all)==1
            fun_predictor = @(s) predictor_initial(var_all,l_all,s,Opt);
        else
            fun_predictor = @(s) predictor_ode(var_all,l_all,s,solver_jacobian,fun);
        end
    else
        error('predictor not set or of unknown type');
    end
    if Opt.predictor_solver
        xi = [var_all(:,end);l_all(end)];
        fun_solve = @(s) merge_arle_pred(fun_predictor,res_arle,s,xi,ds);
        [sp,~,exitflag] = predictor_solver(fun_solve,ds);
        if exitflag<=0
            sp = ds;
        end
    else
        sp = ds;
    end
    xip1 = fun_predictor(sp);
    vp = xip1(1:end-1);
    lp = xip1(end);
end