%% path continuation - stepSize.szyszkowski
%  Adjusts stepsize by the ratio of the angles of the lines connecting the
%  last four consecutive solution points. Also adapts due to needed number 
%  of iterations (see <a href="matlab:doc('stepSize.iterationsPolynomial')">stepSize.iterationsPolynomial</a>).
%  Both adaption factors are weighted by the weights specified in 
%  'weightsSzyszkowski' and then multiplied.
%
%
%   Inputs:
%       oih           -- OptInfoHandle object
%                        
%   Outputs:
%       xi            -- stepsize adaption factor
%
%
%
%  See the <a href="matlab:open('..\doc\html\continuation.html')">documentation</a>. See <a href="matlab:doc('stepSize.control')">other stepsize adaption methods</a>.
%
%   DOI: https://doi.org/10.1007/s004660050513  (adapted version)
%
%   Institute of Dynamics and Vibration Research
%   Leibniz University Hannover
%   17.01.2022 - Tido Kubatschek
%
function [xi] = szyszkowski(oih)
    %% Method of Szyszkowski and Husband
    %
    % collect needed data of Path
    %
    varNeeded = oih.path.varAll(:,end-3:end);
    lNeeded = oih.path.lAll(end-3:end);
    zNeeded = [varNeeded; lNeeded];
    %
    % calculate connecting vectors
    %
    v1 = zNeeded(:,end) - zNeeded(:,end-1);
    v2 = zNeeded(:,end-1) - zNeeded(:,end-2);
    v3 = zNeeded(:,end-2) - zNeeded(:,end-3);
    %
    % calculate angles
    %
    angle = aux.vectorAngle(v1,v2);
    angleM1 = aux.vectorAngle(v2,v3);
    %
    % calculate next angle
    %
    angleP1 = 2*angle - angleM1 * norm(v1)/norm(v2);
    %
    % check if angle is too large
    %
    if angle > oih.opt.stepSizeAngle || angleP1 > oih.opt.stepSizeAngle
        xi = 0.5;
    else
        % correct number of iterations
        %
        if oih.opt.dsMax==inf
            iter = max(oih.solver.output.iterations(end),1);
        else
            iter = oih.solver.output.iterations(end);
        end
        %
        % calculate deviation of iterations
        %
        deviationOfIterations = oih.opt.nIterOpt/iter;
        %
        % get weigths
        weights = oih.opt.weightsSzyszkowski;
        %
        % calculate new step size
        %
        xi = deviationOfIterations^weights(1) *...
            abs(angle/angleP1)^weights(2);
        %
    end
end