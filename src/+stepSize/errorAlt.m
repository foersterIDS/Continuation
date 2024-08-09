%% path continuation - stepSize.errorAlt
%  Adjusts stepsize due to the relative differences of multiple values:
%  -- needed number of iterations and optimal number of iterations
%  -- change of curvature of path and optimal change of curvature (1)
%  -- speed of continuation and optimal speed
%  -- rate of contraction and optimal rate
%  -- distance of predictor to solution point and optimal distance
%  The optimal values are specified by:
%  -- optimal number of iterations: 'nIterOpt'
%  -- optimal speed: 'speedOfContinuation'
%  -- rate of contraction: 'optimalContractionRate'
%  -- distance of predictor: 'predictorDistance'
%  The differences are then devided by the optimal values, weighted by 
%  weights specified in 'weightsMultiplicative' and then added togehter.
%  The sum is then used as an exponent and can be weighted by
%  'stepSizeErrorMax'. Also it is possible to consider the tendency of
%  the adaption in the last step by a PD controller which constants are 
% specified in 'stepSizeErrorPd'.
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
%   20.01.2022 - Tido Kubatschek
%
function [xi] = errorAlt(oih)
    %
    EI = calcError(oih,0);
    %
    %
    K = oih.opt.stepSizeErrorPd;
    %
    if oih.path.nAll > 1 && ~isempty(oih.path.speedOfContinuation) && K(2) > 0
        EIM1 = calcError(oih,1);
    else
        EIM1 = 0;
    end
    %
    E = K(1) * EI + K(2) * (EI - EIM1);
    %
    %% adjustment factor
    %
    xi = 2^(E/oih.opt.stepSizeErrorMax);
end

function EI = calcError(oih,previous)
    %
    %% Weights
    weights = oih.opt.weightsError.';
    %
    %
    % determine length
    lengthPath = oih.path.nAll;
    lengthIterations = length(oih.solver.output.iterations);
    if weights(3) ~= 0
        lengthArrays = length(oih.path.speedOfContinuation);
    elseif weights(4) ~= 0
        lengthArrays = length(oih.solver.output.rateOfConvergence);
    elseif weights(5) ~= 0
        lengthArrays = length(oih.path.xPredictorAll(1,:));
    else
        lengthArrays = 1;
    end
    %
    % determine state
    if previous
        endOfIterations = lengthIterations - 1;
        endOfArray = lengthArrays - 1;
        lengthPath = lengthPath - 1;
    else
        if lengthPath < 3
            endOfIterations = lengthIterations - 1;
            endOfArray = lengthArrays - 1;
        else
            endOfIterations = lengthIterations;
            endOfArray = lengthArrays;
        end
    end
    
    % create vector with oih.path.varAll and oih.path.lAll
    %
    xAll = oih.path.xAll;
    %
    %% define target values
    %
    % wTarget: optimal number of iterations, optimal rate of contraction,
    %           optimal speed of continuation, optimal change of angle,
    %           optimal distance of predictor
    %
    wTarget = [oih.opt.nIterOpt, oih.opt.optimalContractionRate,...
        oih.opt.speedOfContinuation, 1, oih.opt.predictorDistance];
    %
    %
    %% Factor by number of iterations
    % correct number of iterations
    if weights(1) ~= 0 && ~isempty(oih.solver.output.iterations)
        if oih.opt.dsMax==inf
            wIter = max(oih.solver.output.iterations(endOfIterations),1);
        else
            wIter = oih.solver.output.iterations(endOfIterations);
        end
    else
        wIter = wTarget(1);
    end
    %
    %% Factor by contraction rate
    if weights(2) ~= 0 && ~isempty(oih.solver.output.rateOfContraction)
        wContr = oih.solver.output.rateOfContraction(endOfArray);
    else
        wContr = wTarget(2);
    end
    %
    %% Factor by speed of continuation
    if weights(3) ~= 0 && ~isempty(oih.path.speedOfContinuation)
        wSpeed = oih.path.speedOfContinuation(endOfArray);
    else
        wSpeed = wTarget(3);
    end
    %
    %% Factor by change of angle
    %
    % Check if there are enough solution points
    %
    if weights(4) ~= 0 && lengthPath > 3
        %
        % calculate connecting vectors
        %
        v1 = xAll(:,end) - xAll(:,end-1);
        v2 = xAll(:,end-1) - xAll(:,end-2);
        v3 = xAll(:,end-2) - xAll(:,end-3);
        %
        % calculate angles
        %
        angle1 = aux.vectorAngle(v1,v2);
        angle2 = aux.vectorAngle(v2,v3);
        %
        % calculate change of angles
        %
        if angle1 <= oih.opt.solverTol || angle2 <= oih.opt.solverTol
            wCurv = 1;
        else
            wCurv = angle1/angle2;
        end
    else
        wCurv = 1;
    end
    %
    %% Factor by distance of predictor
    %
    if weights(5) ~= 0 && ~isempty(oih.path.xPredictorAll)
        if norm(xAll(:,lengthPath)) > 1e-6
            relDistanceOfPredictor = norm(oih.path.xPredictorAll(:,endOfArray) - xAll(:,lengthPath)) / norm(xAll(:,lengthPath));
            wDist = relDistanceOfPredictor;
        elseif norm(oih.path.xPredictorAll(:,endOfArray)) > 1e-6
            relDistanceOfPredictor = norm(oih.path.xPredictorAll(:,endOfArray) - xAll(:,lengthPath)) / norm(oih.path.xPredictorAll(:,endOfArray));
            wDist = relDistanceOfPredictor;
        else
            wDist = wTarget(5);
        end
    else
        wDist = wTarget(5);
    end
    %
    %% values
    %
    wI = [wIter, wContr, wSpeed, wCurv, wDist];
    %
    %% errors
    %
    eI = (wTarget - wI)./wTarget;    
    %
    %% weighted error sum
    %
    EI = (eI * weights)/sum(weights);
end