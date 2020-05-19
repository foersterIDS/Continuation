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
            options = optimoptions('fsolve','display','off','SpecifyObjectiveGradient',1);
        else
            options = optimoptions('fsolve','display','off');
        end
        solver = @(fun,x0) fsolve(fun,x0,options);
    elseif Opt.solver.fmincon
        %% fmincon
        error('not implemented yet');
        options = optimoptions('fmincon','display','off');
        solver = @(fun,x0) fmincon(fun,x0,[],[],[],[],[],[],[],options);
    elseif Opt.solver.lsqnonlin
        %% lsqnonlin
        if Opt.jacobian
            options = optimoptions('lsqnonlin','display','off','SpecifyObjectiveGradient',1);
        else
            options = optimoptions('lsqnonlin','display','off');
        end
        solver = @(fun,x0) lsqnonlin(fun,x0,[],[],options);
    else
        %% error
        error('No such solver');
    end
end