%% path continuation - corrector.residualSphere
%
%   Institute of Dynamics and Vibration Research
%   Leibniz University Hannover
%   19.05.2020 - Alwin Förster
%
function [residual,jacobian] = residualSphere(x,xAll,ds)
    residual = ds^2-(x-xAll(:,end))'*(x-xAll(:,end));
    jacobian = -2*(x-xAll(:,end))';
end