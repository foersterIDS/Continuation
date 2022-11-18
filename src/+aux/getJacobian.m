%% path continuation - aux.getJacobian
%
%   Institute of Dynamics and Vibration Research
%   Leibniz University Hannover
%   26.05.2020 - Alwin Förster
%
function [jacobian] = getJacobian(fun,v,l,Opt)
    try
        [~,jacobian] = fun(v,l);
    catch
        jacobianV = aux.numericJacobian(@(v) fun(v,l), v, 'diffquot', Opt.diffquot);
        jacobianL = aux.numericJacobian(@(l) fun(v,l), l, 'diffquot', Opt.diffquot);
        jacobian = [jacobianV,jacobianL];
    end
end