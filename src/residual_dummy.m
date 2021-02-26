%% path continuation - residual_dummy
%
%   Institute of Dynamics and Vibration Research
%   Leibniz University Hannover
%   15.02.2021 - Tido Kubatschek
%
function [residual,jacobian] = residual_dummy(x)
    residual = x-1;
    jacobian = 1;
end