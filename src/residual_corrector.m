%% path continuation - residual_corrector
%
%   Institute of Dynamics and Vibration Research
%   Leibniz University Hannover
%   08.05.2020 - Alwin Förster
%   17.09.2020 - Tido Kubatschek
%
function [residual] = residual_corrector(Opt)
    if Opt.corrector.orthogonal
        residual = @(x,x_all,ds,Jac) residual_corrector_orthogonal(x,x_all,ds,Jac,Opt);
    elseif Opt.corrector.sphere
        residual = @(x,x_all,ds) residual_corrector_sphere(x,x_all,ds);
    elseif Opt.corrector.ellipsoid
        RnR = RnRotation([1;0]);
        residual = @(x,x_all,ds) residual_corrector_ellipsoid(x,x_all,ds,RnR);
    elseif Opt.corrector.ellipsoid2
        residual = @(x,x_all,ds) residual_corrector_ellipsoid2(x,x_all,ds);
    elseif Opt.corrector.unique
        residual = @(x,x_all,ds) residual_corrector_unique(x,x_all,ds,Opt);
    elseif Opt.corrector.paraboloid
        residual = @(x,x_all,ds) residual_corrector_paraboloid(x,x_all,ds,Opt);
    else
        error('No such corrector-method');
    end
end