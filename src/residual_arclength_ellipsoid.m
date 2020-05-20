%% path continuation - residual_arclength_sphere
%
%   Institute of Dynamics and Vibration Research
%   Leibniz University Hannover
%   19.05.2020 - Alwin F�rster
%
function [residual,jacobian] = residual_arclength_ellipsoid(x,xs,ds,RnR)
    f = 0.25;%0.05;
    r = ds*[1;f*ones(length(x)-1,1)];
    T = RnR.getTR(diff_xs(xs));
    residual = sum((T*(x-xs(:,end))).^2./r.^2)-1;
    jacobian = (2*(T*(x-xs(:,end)))./r.^2)'*T;
end