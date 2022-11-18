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
%       Solver.output -- contains information of solver, such as the 
%                        needed number of iterations and rate of contraction.
%       Path          -- contains the solution points of the path and the
%                        predictors.
%       Opt           -- contains user inputs, such as optimal values, the
%                        weights, accessible by 'weightsMultiplicative',
%                        the sum weight ('stepSizeErrorMax') and the 
%                        PD constants ('stepSizeErrorPd').
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
function [xi] = errorAlt(Solver,Path,Opt)
    %
    EI = calcError(Solver,Path,Opt,0);
    %
    %
    K = Opt.stepSizeErrorPd;
    %
    if length(Path.lAll) > 1 && ~isempty(Path.speedOfContinuation) && K(2) > 0
        EIM1 = calcError(Solver,Path,Opt,1);
    else
        EIM1 = 0;
    end
    %
    E = K(1) * EI + K(2) * (EI - EIM1);
    %
    %% adjustment factor
    %
    xi = 2^(E/Opt.stepSizeErrorMax);
end

function EI = calcError(Solver,Path,Opt,previous)
    %
    %% Weights
    weights = Opt.weightsError.';
    %
    %
    % determine length
    lengthPath = length(Path.lAll);
    lengthIterations = length(Solver.output.iterations);
    if weights(3) ~= 0
        lengthArrays = length(Path.speedOfContinuation);
    elseif weights(4) ~= 0
        lengthArrays = length(Solver.output.rateOfConvergence);
    elseif weights(5) ~= 0
        lengthArrays = length(Path.xPredictor(1,:));
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
    
    % create vector with Path.varAll and Path.lAll
    %
    xAll = [Path.varAll(:,1:lengthPath);Path.lAll(1:lengthPath)];
    %
    %% define target values
    %
    % wTarget: optimal number of iterations, optimal rate of contraction,
    %           optimal speed of continuation, optimal change of angle,
    %           optimal distance of predictor
    %
    wTarget = [Opt.nIterOpt, Opt.optimalContractionRate,...
        Opt.speedOfContinuation, 1, Opt.predictorDistance];
    %
    %
    %% Factor by number of iterations
    % correct number of iterations
    if weights(1) ~= 0 && ~isempty(Solver.output.iterations)
        if Opt.dsMax==inf
            wIter = max(Solver.output.iterations(endOfIterations),1);
        else
            wIter = Solver.output.iterations(endOfIterations);
        end
    else
        wIter = wTarget(1);
    end
    %
    %% Factor by contraction rate
    if weights(2) ~= 0 && ~isempty(Solver.output.rateOfContraction)
        wContr = Solver.output.rateOfContraction(endOfArray);
    else
        wContr = wTarget(2);
    end
    %
    %% Factor by speed of continuation
    if weights(3) ~= 0 && ~isempty(Path.speedOfContinuation)
        wSpeed = Path.speedOfContinuation(endOfArray);
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
        if angle1 <= Opt.solverTol || angle2 <= Opt.solverTol
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
    if weights(5) ~= 0 && ~isempty(Path.xPredictor)
        if norm(xAll(:,lengthPath)) > 1e-6
            relDistanceOfPredictor = norm(Path.xPredictor(:,endOfArray) - xAll(:,lengthPath)) / norm(xAll(:,lengthPath));
            wDist = relDistanceOfPredictor;
        elseif norm(Path.xPredictor(:,endOfArray)) > 1e-6
            relDistanceOfPredictor = norm(Path.xPredictor(:,endOfArray) - xAll(:,lengthPath)) / norm(Path.xPredictor(:,endOfArray));
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