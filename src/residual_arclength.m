%% path continuation - residual_arclength
%
%   Institute of Dynamics and Vibration Research
%   Leibniz University Hannover
%   08.05.2020 - Alwin Förster
%   17.09.2020 - Tido Kubatschek
%
function [residual] = residual_arclength(Opt)
    if Opt.arclength.linear
        residual = @(x,xs,ds) residual_arclength_linear(x,xs,ds,Opt);
    elseif Opt.arclength.sphere
        residual = @(x,xs,ds) residual_arclength_sphere(x,xs,ds);
    elseif Opt.arclength.ellipsoid
        RnR = RnRotation([1;0]);
        residual = @(x,xs,ds) residual_arclength_ellipsoid(x,xs,ds,RnR);
    elseif Opt.arclength.ellipsoid2
        residual = @(x,xs,ds) residual_arclength_ellipsoid2(x,xs,ds);
    elseif Opt.arclength.unique
        residual = @(x,xs,ds) residual_arclength_unique(x,xs,ds,Opt);
    else
        error('No such arclength-method');
    end
end