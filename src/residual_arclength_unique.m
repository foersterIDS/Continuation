%% path continuation - residual_arclength_unique
%
%   Institute of Dynamics and Vibration Research
%   Leibniz University Hannover
%   22.09.2020 - Alwin Förster
%
function [residual,jacobian] = residual_arclength_unique(x,xs,ds,Opt)
    [a,b] = size(xs);
    residual = x(end)-(xs(end,b)+ds*Opt.direction);
    jacobian = [zeros(1,a-1),1];
end