%% path continuation - residual_fixed_value
%
%   Institute of Dynamics and Vibration Research
%   Leibniz University Hannover
%   20.05.2020 - Alwin Förster
%
function [varargout] = residual_fixed_value(fun,v,l_fix,Opt)
    if Opt.jacobian
        [R,J] = fun(v,l_fix);
        varargout{1} = R;
        varargout{2} = J(1:length(v),1:length(v));
    else
        varargout{1} = fun(v,l_fix);
    end
end