%% path continuation - residual_corrector_sphere
%
%   Institute of Dynamics and Vibration Research
%   Leibniz University Hannover
%   19.05.2020 - Alwin Förster
%
function [residual,jacobian] = residual_corrector_sphere(x,xs,ds)
    residual = ds^2-(x-xs(:,end))'*(x-xs(:,end));
    jacobian = -2*(x-xs(:,end))';
end