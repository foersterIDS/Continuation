%% path continuation - corrector.residualOrthogonal2
%
%   Institute of Dynamics and Vibration Research
%   Leibniz University Hannover
%   14.07.2022 - Alwin Förster
%
function [residual,jacobian] = residualOrthogonal2(x,xAll,ds,fun,Jac,oih)
    [nd,nl] = size(xAll);
    if nl == 1
        % determine method to calculate tangent
        if oih.opt.correctorOrthogonalMethod.secant
            % approximate tangent with secant
            if numel(oih.opt.direction)==1
                tangent = [zeros(nd-1,1);1];
                tangent = oih.opt.direction * ds * tangent;
            else
                tangent = oih.opt.direction * ds;
            end
            xiPredictor = xAll(:,end) + tangent;
        else
            % calc tangent via jacobian
            PathHelp.varAll = xAll(1:end-1,:);
            PathHelp.lAll = xAll(end,:);
            [xiPredictor,tangent] = predictor.initial(PathHelp,ds,oih);
        end
    else
        if oih.opt.correctorOrthogonalMethod.secant
            % approximate tangent with secant
            tangent = xAll(:,end) - xAll(:,end-1);
            xi = xAll(:,end);
            xiPredictor = xi + ds*tangent/norm(tangent);
        else
            % calc tangent via jacobian
            PathHelp.varAll = xAll(1:end-1,:);
            PathHelp.lAll = xAll(end,:);
            [xiPredictor,tangent] = predictor.ode(PathHelp,ds,Jac,[],oih);
        end       
    end
    % jacobian at x
    if oih.opt.jacobian
        [fx,Jx] = fun(x(1:(end-1)),x(end));
    else
        [Jx,fx] = aux.numericJacobian(@(x) fun(x(1:(end-1)),x(end)),x);
    end
    % normal vector to orthogonal hyper plane
    nx = numel(x);
    [~,ind] = min(abs(sqrt(sum(Jx.^2))-median(sqrt(sum(Jx.^2)))));
    indR = 1:nx;
    indR(ind) = [];
    normalRed = Jx(:,indR)\-Jx(:,ind);
    normal = zeros(nx,1);
    normal(ind) = 1;
    normal(indR) = normalRed;
    residual = normal.' * (x - xiPredictor);
    jacobian = normal.'; % Das ist eigentlich falsch, weil Jx von x abhaengt und damit auch normal
end