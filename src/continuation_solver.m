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
            options = optimoptions('fsolve','display','off','SpecifyObjectiveGradient',true,'FunctionTolerance',Opt.solver_tol/10);
        else
            options = optimoptions('fsolve','display','off','SpecifyObjectiveGradient',false,'FunctionTolerance',Opt.solver_tol/10);
        end
        optfun = @(dscale) optimoptions(options,'TypicalX',dscale);
        solver = @(fun,x0,dscale) fsolve(fun,x0,optfun(dscale));
    elseif Opt.solver.lsqnonlin
        %% lsqnonlin
        if Opt.jacobian
            options = optimoptions('lsqnonlin','display','off','SpecifyObjectiveGradient',true,'FunctionTolerance',Opt.solver_tol/10);
        else
            options = optimoptions('lsqnonlin','display','off','SpecifyObjectiveGradient',false,'FunctionTolerance',Opt.solver_tol/10);
        end
        optfun = @(dscale) optimoptions(options,'TypicalX',dscale);
        solver = @(fun,x0,dscale) solver_lsqnonlin(fun,x0,optfun(dscale));
    elseif Opt.solver.newton
        %% basic newton solver
        solver = @(fun,x0,dscale) solver_basic_newton(fun,x0,dscale,Opt);
    else
        %% error
        error('No such solver');
    end
    default_solver_output = struct('iterations',inf,...
                                   'algorithm','non',...
                                   'tolerance',inf);
end