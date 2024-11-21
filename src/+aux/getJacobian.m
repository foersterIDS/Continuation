%% path continuation - aux.getJacobian
%
%   Institute of Dynamics and Vibration Research
%   Leibniz University Hannover
%   26.05.2020 - Alwin Förster
%
function [jacobian] = getJacobian(fun,v,l,oih)
    try
        [~,jacobian] = fun(v,l);
    catch
        aux.printLine(oih,'----> Unable to evaluate user defined jacobian. Using numeric jacobian instead.\n');
        jacobianV = aux.numericJacobian(@(v) fun(v,l), v, 'diffquot', obj.opt.diffquot,'diffStep',oih.opt.diffStep);
        jacobianL = aux.numericJacobian(@(l) fun(v,l), l, 'diffquot', obj.opt.diffquot,'diffStep',oih.opt.diffStep);
        jacobian = [jacobianV,jacobianL];
    end
end