%% path continuation - corrector.residual_unique
%
%   Institute of Dynamics and Vibration Research
%   Leibniz University Hannover
%   22.09.2020 - Alwin Förster
%
function [residual,jacobian] = residual_unique(x,x_all,ds,Opt)
    [a,b] = size(x_all);
    residual = x(end)-(x_all(end,b)+ds*sign(Opt.direction(end)));
    jacobian = [zeros(1,a-1),1];
end