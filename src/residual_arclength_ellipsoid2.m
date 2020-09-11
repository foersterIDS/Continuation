%% path continuation - residual_arclength_ellipsoid2
%
%   Institute of Dynamics and Vibration Research
%   Leibniz University Hannover
%   11.09.2020 - Alwin Förster
%
function [residual,jacobian] = residual_arclength_ellipsoid2(x,xs,ds)
    r = ds;
    if length(xs(1,:))==1
        f = ones(size(xs(:,1)));
    else
        v0 = 10^-15;
        max_mag = max(abs(xs(:,end)-xs(:,end-1))+v0);
        f = (max_mag./(abs(xs(:,end)-xs(:,end-1))+v0));
    end
    residual = sum((f.*(x-xs(:,end))).^2./r.^2)-1;
    jacobian = (2*(f.*(x-xs(:,end)))./r.^2)';
end