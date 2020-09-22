%% path continuation - residual_arclength_ellipsoid2
%
%   Institute of Dynamics and Vibration Research
%   Leibniz University Hannover
%   11.09.2020 - Alwin Förster
%
function [residual,jacobian] = residual_arclength_ellipsoid2(x,xs,ds)
    r = ds;
    fm = 10^6;
    if length(xs(1,:))==1
        f = ones(size(xs(:,1)));
    else
        v0 = 10^-15*ones(size(xs(:,end)));
        axs = abs(xs(:,end));
        adxs = abs(xs(:,end)-xs(:,end-1));
        v0(adxs==0) = axs(adxs==0)./10;
        max_mag = max(adxs+v0);
        nzadxs = adxs+v0;
        f = (max_mag./nzadxs);
        if max(f)~=1
            f = (fm-1)/(max(f)-1)*(f-1)+1;
        end
    end
    residual = sum((f.*(x-xs(:,end))).^2./r.^2)-1;
    jacobian = (2*(f.*(x-xs(:,end)))./r.^2)';
    
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    % So laeuft sich das ganze schnell fest. Gibt es eine bessere Loesung?
    % Kann man das Residuum oder auch die Unbekannten eventuell noch anders
    % skalieren?
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    
end