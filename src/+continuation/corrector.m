%% path continuation - continuation.corrector
%
%   Institute of Dynamics and Vibration Research
%   Leibniz University Hannover
%   08.05.2020 - Alwin Förster
%   17.09.2020 - Tido Kubatschek
%
function [residual] = corrector(fun,Opt)
    if Opt.corrector.orthogonal
        residual = @(x,x_all,ds,Jac) corrector.residual_orthogonal(x,x_all,ds,Jac,Opt);
    elseif Opt.corrector.orthogonal2
        residual = @(x,x_all,ds,Jac) corrector.residual_orthogonal2(x,x_all,ds,fun,Jac,Opt);
    elseif Opt.corrector.sphere
        residual = @(x,x_all,ds,Jac) corrector.residual_sphere(x,x_all,ds);
    elseif Opt.corrector.ellipsoid
        RnR = corrector.RnRotation([1;0]);
        residual = @(x,x_all,ds,Jac) corrector.residual_ellipsoid(x,x_all,ds,RnR);
    elseif Opt.corrector.ellipsoid2
        residual = @(x,x_all,ds,Jac) corrector.residual_ellipsoid2(x,x_all,ds);
    elseif Opt.corrector.unique
        residual = @(x,x_all,ds,Jac) corrector.residual_unique(x,x_all,ds,Opt);
    elseif Opt.corrector.paraboloid
        residual = @(x,x_all,ds,Jac) corrector.residual_paraboloid(x,x_all,ds,Opt);
    else
        error('No such corrector-method');
    end
end