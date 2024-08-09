%% path continuation - corrector.residualOrthogonal
%
%   Institute of Dynamics and Vibration Research
%   Leibniz University Hannover
%   17.09.2020 - Tido Kubatschek
%   22.09.2020 - Alwin Förster
%
function [residual,jacobian] = residualOrthogonal(x,xAll,ds,Jac,oih)
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
    %
    residual = tangent.' * (x - xiPredictor);
    jacobian = tangent.';
end