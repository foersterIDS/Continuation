%% path continuation - corrector.residualEllipsoid
%
%   Institute of Dynamics and Vibration Research
%   Leibniz University Hannover
%   19.05.2020 - Alwin Förster
%
function [residual,jacobian] = residualEllipsoid(x,xAll,ds,RnR)
    f = 0.25;%0.05;
    r = ds*[1;f*ones(length(x)-1,1)];
    T = RnR.getTR(corrector.diffXs(xAll));
    residual = sum((T*(x-xAll(:,end))).^2./r.^2)-1;
    jacobian = (2*(T*(x-xAll(:,end)))./r.^2)'*T;
end