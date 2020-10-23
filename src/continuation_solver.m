%% path continuation - continuation_solver
%
%   Institute of Dynamics and Vibration Research
%   Leibniz University Hannover
%   08.05.2020 - Alwin Förster
%
function [solver,default_solver_output] = continuation_solver(Opt)
    if Opt.solver.fsolve
        %% fsolve
        if Opt.jacobian
            options = optimoptions('fsolve','display','off','SpecifyObjectiveGradient',true,'FunctionTolerance',Opt.solver_tol/2);
        else
            options = optimoptions('fsolve','display','off','SpecifyObjectiveGradient',false,'FunctionTolerance',Opt.solver_tol/2);
        end
        solver = @(fun,x0) fsolve(fun,x0,options);
    elseif Opt.solver.fmincon
        %% fmincon
        warning('fmincon: not implemented yet using fsolve instead.');
%         options = optimoptions('fmincon','display','off');
%         solver = @(fun,x0) fmincon(fun,x0,[],[],[],[],[],[],[],options);
        %% fsolve
        Opt.solve.fsolve = true;
        Opt.solve.fmincon = false;
        if Opt.jacobian
            options = optimoptions('fsolve','display','off','SpecifyObjectiveGradient',true,'FunctionTolerance',Opt.solver_tol/2);
        else
            options = optimoptions('fsolve','display','off','SpecifyObjectiveGradient',false,'FunctionTolerance',Opt.solver_tol/2);
        end
        solver = @(fun,x0) fsolve(fun,x0,options);
    elseif Opt.solver.lsqnonlin
        %% lsqnonlin
        if Opt.jacobian
            options = optimoptions('lsqnonlin','display','off','SpecifyObjectiveGradient',true,'FunctionTolerance',Opt.solver_tol,'OptimalityTolerance',Opt.solver_tol);
        else
            options = optimoptions('lsqnonlin','display','off','SpecifyObjectiveGradient',false,'FunctionTolerance',Opt.solver_tol,'OptimalityTolerance',Opt.solver_tol);
        end
        solver = @(fun,x0) lsqnonlin(fun,x0,[],[],options);
    elseif Opt.solver.newton
        %% basic newton solver
        solver = @(fun,x0) basic_newton_solver(fun,x0,Opt);
    else
        %% error
        error('No such solver');
    end
    default_solver_output = struct('iterations',inf,...
                                   'algorithm','non',...
                                   'tolerance',inf);
end