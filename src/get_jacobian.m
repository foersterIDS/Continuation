%% path continuation - get_jacobian
%
%   Institute of Dynamics and Vibration Research
%   Leibniz University Hannover
%   26.05.2020 - Alwin Förster
%
function [jacobian] = get_jacobian(fun,v,l)
    try
        [~,jacobian] = fun(v,l);
    catch
        jacobian = numeric_jacobian(f, x);
    end
end