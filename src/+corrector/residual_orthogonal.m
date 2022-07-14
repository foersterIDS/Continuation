%% path continuation - corrector.residual_orthogonal
%
%   Institute of Dynamics and Vibration Research
%   Leibniz University Hannover
%   17.09.2020 - Tido Kubatschek
%   22.09.2020 - Alwin Förster
%
function [residual,jacobian] = residual_orthogonal(x,x_all,ds,Jac,Opt)
    [nd,nl] = size(x_all);
    if nl == 1
        % determine method to calculate tangent
        if Opt.corrector_orthogonal_method.secant
            % approximate tangent with secant
            if numel(Opt.direction)==1
                tangent = [zeros(nd-1,1);1];
                tangent = Opt.direction * ds * tangent;
            else
                tangent = Opt.direction * ds;
            end
            xi_predictor = x_all(:,end) + tangent;
        else
            % calc tangent via jacobian
            Path_help.var_all = x_all(1:end-1,:);
            Path_help.l_all = x_all(end,:);
            [xi_predictor,tangent] = predictor.initial(Path_help,ds,Opt);
        end
    else
        if Opt.corrector_orthogonal_method.secant
            % approximate tangent with secant
            tangent = x_all(:,end) - x_all(:,end-1);
            xi = x_all(:,end);
            xi_predictor = xi + ds*tangent/norm(tangent);
        else
            % calc tangent via jacobian
            Path_help.var_all = x_all(1:end-1,:);
            Path_help.l_all = x_all(end,:);
            [xi_predictor,tangent] = predictor.ode(Path_help,ds,Jac,[],Opt);
        end       
    end
    %
    residual = tangent.' * (x - xi_predictor);
    jacobian = tangent.';
end