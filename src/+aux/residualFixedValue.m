%% path continuation - aux.residualFixedValue
%
%   Institute of Dynamics and Vibration Research
%   Leibniz University Hannover
%   20.05.2020 - Alwin Förster
%
function [varargout] = residualFixedValue(func,v,lFix,oih)
    if oih.opt.jacobian
        [R,J] = func(v,lFix);
        [n1,n2] = size(J);
        if n2<=n1 && n2<numel(v)
            Jl = aux.numericJacobian(@(v) func(v,lFix), v, 'derivativeDimensions', (n2+1):numel(v), 'diffquot', oih.opt.diffquot, 'centralValue', R);
        else
            Jl = [];
        end
        J = [J,Jl];
        varargout{1} = R;
        varargout{2} = J(1:numel(v),1:numel(v));
    else
        varargout{1} = func(v,lFix);
    end
end