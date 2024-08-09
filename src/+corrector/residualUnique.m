%% path continuation - corrector.residualUnique
%
%   Institute of Dynamics and Vibration Research
%   Leibniz University Hannover
%   22.09.2020 - Alwin Förster
%
function [residual,jacobian] = residualUnique(x,xAll,ds,oih)
    [a,b] = size(xAll);
    residual = x(end)-(xAll(end,b)+ds*sign(oih.opt.direction(end)));
    jacobian = [zeros(1,a-1),1];
end