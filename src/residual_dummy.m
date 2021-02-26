%% path continuation - residual_dummy
%
%   Institute of Dynamics and Vibration Research
%   Leibniz University Hannover
%   15.02.2021 - Tido Kubatschek
%
function [residual,jacobian] = residual_dummy(x,xs,ds)
    residual = x(end)-1;
    jacobian = zeros(1, numel(x));
    jacobian(end) = 1;
end