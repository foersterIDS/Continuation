%% path continuation - aux.residual_suspend_continuation
%
%   Institute of Dynamics and Vibration Research
%   Leibniz University Hannover
%   20.12.2021 - Alwin Förster
%
function [varargout] = residual_suspend_continuation(func,v,l_fix,Opt)
    nv = numel(v);
    if Opt.jacobian
        [R,J] = func(v,l_fix);
        [n1,n2] = size(J);
        if n2<=n1 && n2<nv
            Jl = numeric_jacobian(@(v) func(v,l_fix), v, 'derivative_dimensions', (n2+1):numel(v), 'diffquot', Opt.diffquot, 'central_value', R);
        else
            J = J(:,1:nv);
            Jl = [];
        end
        J = [J,Jl];
        varargout{1} = R;
        varargout{2} = J;
    else
        varargout{1} = func(v,l_fix);
    end
end