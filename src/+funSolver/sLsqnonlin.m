%% path continuation - fun_solver.s_lsqnonlin
%
%   Institute of Dynamics and Vibration Research
%   Leibniz University Hannover
%   23.10.2020 - Alwin Förster
%
function [x,fval,exitflag,output,jacobian] = s_lsqnonlin(fun,x0,options)
    [x,~,fval,exitflag,output,~,jacobian] = lsqnonlin(fun,x0,[],[],options);
end