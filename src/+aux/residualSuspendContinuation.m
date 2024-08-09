%% path continuation - aux.residualSuspendContinuation
%
%   Institute of Dynamics and Vibration Research
%   Leibniz University Hannover
%   20.12.2021 - Alwin Förster
%
function [varargout] = residualSuspendContinuation(func,v,lFix,oih)
    nv = numel(v);
    if oih.opt.jacobian
        [R,J] = func(v,lFix);
        [n1,n2] = size(J);
        if n2<=n1 && n2<nv
            Jl = numericJacobian(@(v) func(v,lFix), v, 'derivativeDimensions', (n2+1):numel(v), 'diffquot', oih.opt.diffquot, 'centralValue', R);
        else
            J = J(:,1:nv);
            Jl = [];
        end
        J = [J,Jl];
        varargout{1} = R;
        varargout{2} = J;
    else
        varargout{1} = func(v,lFix);
    end
end