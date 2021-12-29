%% path continuation - residual_corrector_ellipsoid2
%
%   Institute of Dynamics and Vibration Research
%   Leibniz University Hannover
%   11.09.2020 - Alwin Förster
%
function [residual,jacobian] = residual_corrector_ellipsoid2(x,x_all,ds)
    r = ds;
    fm = 10^6;
    if length(x_all(1,:))==1
        f = ones(size(x_all(:,1)));
    else
        v0 = 10^-15*ones(size(x_all(:,end)));
        axs = abs(x_all(:,end));
        adxs = abs(x_all(:,end)-x_all(:,end-1));
        v0(adxs==0) = axs(adxs==0)./10;
        max_mag = max(adxs+v0);
        nzadxs = adxs+v0;
        f = (max_mag./nzadxs);
        if max(f)~=1
            f = (fm-1)/(max(f)-1)*(f-1)+1;
        end
    end
    residual = sum((f.*(x-x_all(:,end))).^2./r.^2)-1;
    jacobian = (2*(f.*(x-x_all(:,end)))./r.^2)';
    
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    % So laeuft sich das ganze schnell fest. Gibt es eine bessere Loesung?
    % Kann man das Residuum oder auch die Unbekannten eventuell noch anders
    % skalieren?
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    
end