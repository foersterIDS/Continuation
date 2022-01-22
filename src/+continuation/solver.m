%% path continuation - continuation.solver
%
%   Institute of Dynamics and Vibration Research
%   Leibniz University Hannover
%   08.05.2020 - Alwin Förster
%
function [solver,predictor_solver,num_jac_solver,default_solver_output] = solver(Opt,output_flag)
    predictor_solver = [];
    if output_flag
        outfun = @output_fun;
    else
        outfun = @do_nothing;
    end
    if Opt.solver.fsolve
        %% fsolve
        if Opt.jacobian
            options = optimoptions('fsolve','display','off','SpecifyObjectiveGradient',true,'FunctionTolerance',Opt.solver_tol,'StepTolerance',Opt.solver_tol,'OptimalityTolerance',Opt.solver_tol,'MaxIterations',Opt.solver_max_iterations,'OutputFcn', outfun);
        else
            options = optimoptions('fsolve','display','off','SpecifyObjectiveGradient',false,'FunctionTolerance',Opt.solver_tol,'StepTolerance',Opt.solver_tol,'OptimalityTolerance',Opt.solver_tol,'MaxIterations',Opt.solver_max_iterations,'OutputFcn', outfun);
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
            options = optimoptions('lsqnonlin','display','off','SpecifyObjectiveGradient',true,'FunctionTolerance',Opt.solver_tol,'MaxIterations',Opt.solver_max_iterations,'OutputFcn', outfun);
        else
            options = optimoptions('lsqnonlin','display','off','SpecifyObjectiveGradient',false,'FunctionTolerance',Opt.solver_tol,'MaxIterations',Opt.solver_max_iterations,'OutputFcn', outfun);
        end
        optfun = @(dscale) optimoptions(options,'TypicalX',dscale);
        solver = @(fun,x0,dscale) fun_solver.s_lsqnonlin(fun,x0,optfun(dscale));
        if Opt.predictor_solver
            predictor_options = optimoptions('lsqnonlin','display','off','SpecifyObjectiveGradient',true);
            predictor_solver = @(fun_predictor,s0) fun_solver.s_lsqnonline(fun_predictor,s0,predictor_options);
        end
        opt_num_jac = optimoptions(options,'SpecifyObjectiveGradient',false);
        num_jac_solver = @(fun,x0) fun_solver.s_lsqnonline(fun,x0,opt_num_jac);
    elseif Opt.solver.newton
        %% basic newton solver
        solver = @(fun,x0,dscale) fun_solver.basic_newton(fun,x0,dscale,output_flag,Opt);
        if Opt.predictor_solver
            predictor_Opt = Opt;
            predictor_Opt.jacobian = true;
            predictor_solver = @(fun_predictor,s0) fun_solver.basic_newton(fun_predictor,s0,1,predictor_Opt);
        end
        opt_num_jac = Opt;
        opt_num_jac.jacobian = false;
        num_jac_solver = @(fun,x0) fun_solver.basic_newton(fun,x0,1,opt_num_jac);
    else
        %% error
        error('No such solver');
    end
    default_solver_output = struct('iterations',inf,...
                                   'algorithm','non',...
                                   'tolerance',inf);
end

function stop = output_fun(x,optimValues,state)
    stop = false;
    global solver_stepsizes;
    switch state
        case 'init'
            solver_stepsizes = [];
        case 'iter'
            iter = optimValues.iteration;               % Iteration
            normd = optimValues.stepsize;               % Norm of step
            solver_stepsizes = [solver_stepsizes; iter normd];
    end
end

function stop = do_nothing(x,optimValues,state)
    stop = false;
end
