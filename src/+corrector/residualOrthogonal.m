%% path continuation - corrector.residualOrthogonal
%
%   Institute of Dynamics and Vibration Research
%   Leibniz University Hannover
%   17.09.2020 - Tido Kubatschek
%   22.09.2020 - Alwin Förster
%
function [residual,jacobian] = residualOrthogonal(x,xAll,ds,Jac,Opt)
    [nd,nl] = size(xAll);
    if nl == 1
        % determine method to calculate tangent
        if Opt.correctorOrthogonalMethod.secant
            % approximate tangent with secant
            if numel(Opt.direction)==1
                tangent = [zeros(nd-1,1);1];
                tangent = Opt.direction * ds * tangent;
            else
                tangent = Opt.direction * ds;
            end
            xiPredictor = xAll(:,end) + tangent;
        else
            % calc tangent via jacobian
            PathHelp.varAll = xAll(1:end-1,:);
            PathHelp.lAll = xAll(end,:);
            [xiPredictor,tangent] = predictor.initial(PathHelp,ds,Opt);
        end
    else
        if Opt.correctorOrthogonalMethod.secant
            % approximate tangent with secant
            tangent = xAll(:,end) - xAll(:,end-1);
            xi = xAll(:,end);
            xiPredictor = xi + ds*tangent/norm(tangent);
        else
            % calc tangent via jacobian
            PathHelp.varAll = xAll(1:end-1,:);
            PathHelp.lAll = xAll(end,:);
            [xiPredictor,tangent] = predictor.ode(PathHelp,ds,Jac,[],Opt);
        end       
    end
    %
    residual = tangent.' * (x - xiPredictor);
    jacobian = tangent.';
end