%% path continuation - get_jacobian
%
%   Institute of Dynamics and Vibration Research
%   Leibniz University Hannover
%   26.05.2020 - Alwin Förster
%
function [jacobian] = get_jacobian(fun,v,l,Opt)
    try
        [~,jacobian] = fun(v,l);
    catch
        jacobian_v = numeric_jacobian(@(v) fun(v,l), v, 'diffquot', Opt.diffquot);
        jacobian_l = numeric_jacobian(@(l) fun(v,l), l, 'diffquot', Opt.diffquot);
        jacobian = [jacobian_v,jacobian_l];
    end
end