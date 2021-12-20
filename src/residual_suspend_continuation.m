%% path continuation - residual_suspend_continuation
%
%   Institute of Dynamics and Vibration Research
%   Leibniz University Hannover
%   20.05.2020 - Alwin Förster
%
function [varargout] = residual_suspend_continuation(fun,x,Opt)
    v = x(1:(end-1));
    l_fix = x(end);
    if Opt.jacobian
        [R,J] = fun(v,l_fix);
        [n1,n2] = size(J);
        if n2<=n1 && n2<numel(v)
            Jl = numeric_jacobian(@(v) fun(v,l_fix), v, 'derivative_dimensions', (n2+1):numel(v), 'diffquot', Opt.diffquot, 'central_value', R);
        else
            Jl = [];
        end
        J2 = zeros(size(R));
        R = [R;0];
        J = [J,Jl,J2];
        J = [J;zeros(1,numel(J(1,:)))];
        varargout{1} = R;
        varargout{2} = J;
    else
        varargout{1} = fun(v,l_fix);
    end
end