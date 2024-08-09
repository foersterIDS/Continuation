%% path continuation - stepSize.angleChange
%  Adjusts stepsize by the change of the angles of the lines connecting the
%  last four consecutive solution points. Also adapts due to needed number 
%  of iterations (see <a href="matlab:doc('stepSize.iterationsPolynomial')">stepSize.iterationsPolynomial</a>).
%  Both adaption factors are weighted by the weights specified in 
%  'weightsAngleChange' and then multiplied.
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
%   Institute of Dynamics and Vibration Research
%   Leibniz University Hannover
%   17.01.2022 - Tido Kubatschek
%
function [xi] = angleChange(oih)
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
    angle1 = aux.vectorAngle(v1,v2);
    angle2 = aux.vectorAngle(v2,v3);
    %
    % calculate change of curvature
    %
    changeOfAngle = angle2/angle1;
    %
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
    %
    weights = oih.opt.weightsAngleChange;
    %
    % calculate adaption factor
    %
    xi = deviationOfIterations^weights(1) * changeOfAngle^weights(2);
    %
end