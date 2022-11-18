%% path continuation - aux.mergeResiduals
%
%   Institute of Dynamics and Vibration Research
%   Leibniz University Hannover
%   19.05.2020 - Alwin Förster
%   16.09.2020 - Tido Kubatschek 
%
function [varargout] = mergeResiduals(func,resCorr,x,xAll,ds,lastJacobian,Opt)
    if Opt.jacobian
        [R1,J1] = func(x(1:end-1),x(end));
        [R2,J2] = resCorr(x,xAll,ds,lastJacobian);
        [n11,n12] = size(J1);
        if n12<=n11
            J1l = aux.numericJacobian(@(x) func(x(1:end-1),x(end)), x, 'derivativeDimensions', (n12+1):numel(x), 'diffquot', Opt.diffquot, 'centralValue', R1);
        else
            J1l = [];
        end
        varargout{1} = [R1;R2];
        varargout{2} = [J1,J1l;J2];
    else
        R1 = func(x(1:end-1),x(end));
        R2 = resCorr(x,xAll,ds,lastJacobian);
        varargout{1} = [R1;R2];
    end
end