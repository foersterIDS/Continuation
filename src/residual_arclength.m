%% path continuation - residual_arclength
%
%   Institute of Dynamics and Vibration Research
%   Leibniz University Hannover
%   08.05.2020 - Alwin Förster
%
function [residual] = residual_arclength(Opt)
    if Opt.arclength.linear
        error('Muss noch implementiert werden!');
    elseif Opt.arclength.sphere
        residual = @(x,xs,ds) ds^2-(x-xs(:,end))'*(x-xs(:,end));
    elseif Opt.arclength.ellipsoid
        RnR = RnRotation([1;0]);
        f = 0.25;%0.05;
        r = @(x,ds) ds*[1;f*ones(length(x)-1,1)];
        residual = @(x,xs,ds) sum((RnR.getTR(diff_xs(xs))*(x-xs(:,end))).^2./r(x,ds).^2)-1;
    else
        error('No such arclength-method');
    end
end