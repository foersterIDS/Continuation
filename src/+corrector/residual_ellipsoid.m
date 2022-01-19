%% path continuation - corrector.residual_ellipsoid
%
%   Institute of Dynamics and Vibration Research
%   Leibniz University Hannover
%   19.05.2020 - Alwin Förster
%
function [residual,jacobian] = residual_ellipsoid(x,x_all,ds,RnR)
    f = 0.25;%0.05;
    r = ds*[1;f*ones(length(x)-1,1)];
    T = RnR.getTR(corrector.diff_xs(x_all));
    residual = sum((T*(x-x_all(:,end))).^2./r.^2)-1;
    jacobian = (2*(T*(x-x_all(:,end)))./r.^2)'*T;
end