%% path continuation - continuation.corrector
%
%   Institute of Dynamics and Vibration Research
%   Leibniz University Hannover
%   08.05.2020 - Alwin Förster
%   17.09.2020 - Tido Kubatschek
%
function [residual] = corrector(fun,Opt)
    if Opt.corrector.orthogonal
        residual = @(x,xAll,ds,Jac) corrector.residualOrthogonal(x,xAll,ds,Jac,Opt);
    elseif Opt.corrector.orthogonal2
        residual = @(x,xAll,ds,Jac) corrector.residualOrthogonal2(x,xAll,ds,fun,Jac,Opt);
    elseif Opt.corrector.sphere
        residual = @(x,xAll,ds,Jac) corrector.residualSphere(x,xAll,ds);
    elseif Opt.corrector.ellipsoid
        RnR = corrector.RnRotation([1;0]);
        residual = @(x,xAll,ds,Jac) corrector.residualEllipsoid(x,xAll,ds,RnR);
    elseif Opt.corrector.ellipsoid2
        residual = @(x,xAll,ds,Jac) corrector.residualEllipsoid2(x,xAll,ds);
    elseif Opt.corrector.unique
        residual = @(x,xAll,ds,Jac) corrector.residualUnique(x,xAll,ds,Opt);
    elseif Opt.corrector.paraboloid
        residual = @(x,xAll,ds,Jac) corrector.residualParaboloid(x,xAll,ds,Opt);
    else
        error('No such corrector-method');
    end
end