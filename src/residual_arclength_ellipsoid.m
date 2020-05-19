%% path continuation - residual_arclength_sphere
%
%   Institute of Dynamics and Vibration Research
%   Leibniz University Hannover
%   19.05.2020 - Alwin Förster
%
function [residual] = residual_arclength_ellipsoid(x,xs,ds)
    RnR = RnRotation([1;0]);
    f = 0.25;%0.05;
    r = @(x,ds) ds*[1;f*ones(length(x)-1,1)];
    residual = sum((RnR.getTR(diff_xs(xs))*(x-xs(:,end))).^2./r(x,ds).^2)-1;
    % TODO: jacobian
end