%% path continuation - corrector.residualUnique
%
%   Institute of Dynamics and Vibration Research
%   Leibniz University Hannover
%   22.09.2020 - Alwin Förster
%
function [residual,jacobian] = residualUnique(x,xAll,ds,Opt)
    [a,b] = size(xAll);
    residual = x(end)-(xAll(end,b)+ds*sign(Opt.direction(end)));
    jacobian = [zeros(1,a-1),1];
end