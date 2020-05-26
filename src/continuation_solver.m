%% path continuation - continuation_solver
%
%   Institute of Dynamics and Vibration Research
%   Leibniz University Hannover
%   08.05.2020 - Alwin Förster
%
function [solver] = continuation_solver(Opt)
    if Opt.solver.fsolve
        %% fsolve
        if Opt.jacobian
            options = optimoptions('fsolve','display','off','SpecifyObjectiveGradient',true);
        else
            options = optimoptions('fsolve','display','off','SpecifyObjectiveGradient',false);
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
            options = optimoptions('fsolve','display','off','SpecifyObjectiveGradient',true);
        else
            options = optimoptions('fsolve','display','off','SpecifyObjectiveGradient',false);
        end
        solver = @(fun,x0) fsolve(fun,x0,options);
    elseif Opt.solver.lsqnonlin
        %% lsqnonlin
        if Opt.jacobian
            options = optimoptions('lsqnonlin','display','off','SpecifyObjectiveGradient',true);
        else
            options = optimoptions('lsqnonlin','display','off','SpecifyObjectiveGradient',false);
        end
        solver = @(fun,x0) lsqnonlin(fun,x0,[],[],options);
    elseif Opt.solver.newton
        %% basic newton solver
        solver = @(fun,x0) basic_newton_solver(fun,x0,Opt);
    else
        %% error
        error('No such solver');
    end
end