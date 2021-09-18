%% path continuation - continuation_solver
%
%   Institute of Dynamics and Vibration Research
%   Leibniz University Hannover
%   08.05.2020 - Alwin Förster
%
function [solver,predictor_solver,num_jac_solver,default_solver_output] = continuation_solver(Opt)
    predictor_solver = [];
    if Opt.solver.fsolve
        %% fsolve
        if Opt.jacobian
            options = optimoptions('fsolve','display','off','SpecifyObjectiveGradient',true,'FunctionTolerance',Opt.solver_tol/10,'StepTolerance',Opt.solver_tol/10,'OptimalityTolerance',Opt.solver_tol/10);
        else
            options = optimoptions('fsolve','display','off','SpecifyObjectiveGradient',false,'FunctionTolerance',Opt.solver_tol/10,'StepTolerance',Opt.solver_tol/10,'OptimalityTolerance',Opt.solver_tol/10);
        end
        optfun = @(dscale) optimoptions(options,'TypicalX',dscale);
        solver = @(fun,x0,dscale) fsolve(fun,x0,optfun(dscale));
        if Opt.predictor_solver
            predictor_options = optimoptions('fsolve','display','off','SpecifyObjectiveGradient',true);
            predictor_solver = @(fun_predictor,s0) fsolve(fun_predictor,s0,predictor_options);
        end
        opt_num_jac = optimoptions(options,'SpecifyObjectiveGradient',false);
        num_jac_solver = @(fun,x0) fsolve(fun,x0,opt_num_jac);
    elseif Opt.solver.lsqnonlin
        %% lsqnonlin
        if Opt.jacobian
            options = optimoptions('lsqnonlin','display','off','SpecifyObjectiveGradient',true,'FunctionTolerance',Opt.solver_tol/10);
        else
            options = optimoptions('lsqnonlin','display','off','SpecifyObjectiveGradient',false,'FunctionTolerance',Opt.solver_tol/10);
        end
        optfun = @(dscale) optimoptions(options,'TypicalX',dscale);
        solver = @(fun,x0,dscale) solver_lsqnonlin(fun,x0,optfun(dscale));
        if Opt.predictor_solver
            predictor_options = optimoptions('lsqnonlin','display','off','SpecifyObjectiveGradient',true);
            predictor_solver = @(fun_predictor,s0) solver_lsqnonlin(fun_predictor,s0,predictor_options);
        end
        opt_num_jac = optimoptions(options,'SpecifyObjectiveGradient',false);
        num_jac_solver = @(fun,x0) solver_lsqnonlin(fun,x0,opt_num_jac);
    elseif Opt.solver.newton
        %% basic newton solver
        solver = @(fun,x0,dscale) solver_basic_newton(fun,x0,dscale,Opt);
        if Opt.predictor_solver
            predictor_Opt = Opt;
            predictor_Opt.jacobian = true;
            predictor_solver = @(fun_predictor,s0) solver_basic_newton(fun_predictor,s0,1,predictor_Opt);
        end
        opt_num_jac = Opt;
        opt_num_jac.jacobian = false;
        num_jac_solver = @(fun,x0) solver_basic_newton(fun,x0,1,opt_num_jac);
    else
        %% error
        error('No such solver');
    end
    default_solver_output = struct('iterations',inf,...
                                   'algorithm','non',...
                                   'tolerance',inf);
end