%% path continuation - corrector.residual_orthogonal2
%
%   Institute of Dynamics and Vibration Research
%   Leibniz University Hannover
%   14.07.2022 - Alwin Förster
%
function [residual,jacobian] = residual_orthogonal2(x,x_all,ds,fun,Jac,Opt)
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
    % jacobian at x
    if Opt.jacobian
        [fx,Jx] = fun(x(1:(end-1)),x(end));
    else
        [Jx,fx] = aux.numeric_jacobian(@(x) fun(x(1:(end-1)),x(end)),x);
    end
    % normal vector to orthogonal hyper plane
    nx = numel(x);
    [~,ind] = min(abs(sqrt(sum(Jx.^2))-median(sqrt(sum(Jx.^2)))));
    ind_r = 1:nx;
    ind_r(ind) = [];
    normal_red = Jx(:,ind_r)\-Jx(:,ind);
    normal = zeros(nx,1);
    normal(ind) = 1;
    normal(ind_r) = normal_red;
    residual = normal.' * (x - xi_predictor);
    jacobian = normal.'; % Das ist eigentlich falsch, weil Jx von x abhaengt und damit auch normal
end