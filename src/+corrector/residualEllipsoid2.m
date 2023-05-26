%% path continuation - corrector.residualEllipsoid2
%
%   Institute of Dynamics and Vibration Research
%   Leibniz University Hannover
%   11.09.2020 - Alwin Förster
%
function [residual,jacobian] = residualEllipsoid2(x,xAll,ds)
    r = ds;
    fm = 10^6;
    if length(xAll(1,:))==1
        f = ones(size(xAll(:,1)));
    else
        v0 = 10^-15*ones(size(xAll(:,end)));
        axs = abs(xAll(:,end));
        adxs = abs(xAll(:,end)-xAll(:,end-1));
        v0(adxs==0) = axs(adxs==0)./10;
        maxMag = max(adxs+v0);
        nzadxs = adxs+v0;
        f = (maxMag./nzadxs);
        if max(f)~=1
            f = (fm-1)/(max(f)-1)*(f-1)+1;
        end
    end
    residual = sum((f.*(x-xAll(:,end))).^2./r.^2)-1;
    jacobian = (2*(f.*(x-xAll(:,end)))./r.^2)';
    
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    % So laeuft sich das ganze schnell fest. Gibt es eine bessere Loesung?
    % Kann man das Residuum oder auch die Unbekannten eventuell noch anders
    % skalieren?
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    
end