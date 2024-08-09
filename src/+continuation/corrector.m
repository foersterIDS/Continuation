%% path continuation - continuation.corrector
%
%   Institute of Dynamics and Vibration Research
%   Leibniz University Hannover
%   08.05.2020 - Alwin Förster
%   17.09.2020 - Tido Kubatschek
%
function [residual] = corrector(fun,oih)
    if oih.opt.corrector.orthogonal
        residual = @(x,xAll,ds,Jac) corrector.residualOrthogonal(x,xAll,ds,Jac,oih);
    elseif oih.opt.corrector.orthogonal2
        residual = @(x,xAll,ds,Jac) corrector.residualOrthogonal2(x,xAll,ds,fun,Jac,oih);
    elseif oih.opt.corrector.sphere
        residual = @(x,xAll,ds,Jac) corrector.residualSphere(x,xAll,ds);
    elseif oih.opt.corrector.ellipsoid
        RnR = corrector.RnRotation([1;0]);
        residual = @(x,xAll,ds,Jac) corrector.residualEllipsoid(x,xAll,ds,RnR);
    elseif oih.opt.corrector.ellipsoid2
        residual = @(x,xAll,ds,Jac) corrector.residualEllipsoid2(x,xAll,ds);
    elseif oih.opt.corrector.unique
        residual = @(x,xAll,ds,Jac) corrector.residualUnique(x,xAll,ds,oih);
    elseif oih.opt.corrector.paraboloid
        residual = @(x,xAll,ds,Jac) corrector.residualParaboloid(x,xAll,ds,oih);
    else
        error('No such corrector-method');
    end
end