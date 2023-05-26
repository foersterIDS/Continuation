%% path continuation - funSolver.sLsqnonlin
%
%   Institute of Dynamics and Vibration Research
%   Leibniz University Hannover
%   23.10.2020 - Alwin Förster
%
function [x,fval,exitflag,output,jacobian] = sLsqnonlin(fun,x0,options)
    [x,~,fval,exitflag,output,~,jacobian] = lsqnonlin(fun,x0,[],[],options);
end