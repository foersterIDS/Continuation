%% path continuation - residual_corrector_sphere
%
%   Institute of Dynamics and Vibration Research
%   Leibniz University Hannover
%   19.05.2020 - Alwin Förster
%
function [residual,jacobian] = residual_corrector_sphere(x,x_all,ds)
    residual = ds^2-(x-x_all(:,end))'*(x-x_all(:,end));
    jacobian = -2*(x-x_all(:,end))';
end