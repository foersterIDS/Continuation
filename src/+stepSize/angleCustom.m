%% path continuation - stepSize.angleCustom
%  Adjusts stepsize by the ratio of the angle of the lines connecting 
%  three consecutive solution points and an optimal angle specified in 
%  'stepSizeAngle'. Also adapts due to needed number of iterations 
%  (see <a href="matlab:doc('stepSize.iterationsPolynomial')">stepSize.iterationsPolynomial</a>).
%  Both adaption factors are weighted by the weights specified in 
%  'weightsAngleCustom' and then multiplied.
%
%
%   Inputs:
%       Solver.output -- contains information of solver, such as the 
%                        needed number of iterations.
%       Path          -- contains the solution points of the path
%       Opt           -- contains user inputs, such as optimal number of
%                        iterations, accessible by 'nIterOpt' and the
%                        weights, accessible by 'weightsAngleCustom'.
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
function [xi] = angleCustom(Solver,Path,Opt)
    %
    % collect needed data of Path
    %
    varNeeded = Path.varAll(:,end-2:end);
    lNeeded = Path.lAll(end-2:end);
    zNeeded = [varNeeded;lNeeded];
    %
    % calculate connecting vectors
    %
    v1 = zNeeded(:,end) - zNeeded(:,end-1);
    v2 = zNeeded(:,end-1) - zNeeded(:,end-2);
    %
    % calculate angle
    %
    angle = aux.vectorAngle(v1,v2);
    %
    % correct number of iterations
    %
    if Opt.dsMax==inf
        iter = max(Solver.output.iterations(end),1);
    else
        iter = Solver.output.iterations(end);
    end
    %
    % calculate deviation of iterations
    %
    deviationOfIterations = Opt.nIterOpt/iter;
    %
    % calculate ratio
    %
    ratio = angle / Opt.stepSizeAngle;
    %
    if ratio < 0.9 % if ratio < 0.9, set ratio to 0.9
        ratio = 0.9;
    elseif ratio > 1.2 % high ratios are punished stronger (but ratioMax = 5)
        ratio = 5;
    end
    %
    % get weigths
    weights = Opt.weightsAngleCustom;
    %
    xi = deviationOfIterations^weights(1) * (1/ratio)^weights(2);
    %
end