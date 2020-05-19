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
        residual = @(x,xs,ds) residual_arclength_sphere(x,xs,ds);
    elseif Opt.arclength.ellipsoid
        residual = @(x,xs,ds) residual_arclength_ellipsoid(x,xs,ds);
        Opt.jacobian = false;
        warning('Jacobian is not used due to ellipsoid');
    else
        error('No such arclength-method');
    end
end