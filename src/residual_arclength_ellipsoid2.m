%% path continuation - residual_arclength_ellipsoid2
%
%   Institute of Dynamics and Vibration Research
%   Leibniz University Hannover
%   11.09.2020 - Alwin Förster
%
function [residual,jacobian] = residual_arclength_ellipsoid2(x,xs,ds)
    max_mag = max(abs(xs(:,end))+10^-15);
    r = ds*(abs(xs(:,end))/max_mag);
    residual = sum(((x-xs(:,end))).^2./r.^2)-1;
    jacobian = (2*((x-xs(:,end)))./r.^2)';
end