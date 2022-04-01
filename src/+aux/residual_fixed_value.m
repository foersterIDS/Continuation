%% path continuation - aux.residual_fixed_value
%
%   Institute of Dynamics and Vibration Research
%   Leibniz University Hannover
%   20.05.2020 - Alwin Förster
%
function [varargout] = residual_fixed_value(func,v,l_fix,Opt)
    if Opt.jacobian
        [R,J] = func(v,l_fix);
        [n1,n2] = size(J);
        if n2<=n1 && n2<numel(v)
            Jl = aux.numeric_jacobian(@(v) func(v,l_fix), v, 'derivative_dimensions', (n2+1):numel(v), 'diffquot', Opt.diffquot, 'central_value', R);
        else
            Jl = [];
        end
        J = [J,Jl];
        varargout{1} = R;
        varargout{2} = J(1:numel(v),1:numel(v));
    else
        varargout{1} = func(v,l_fix);
    end
end