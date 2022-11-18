%% path continuation - stepSize.multiplicativeAlt
%  Adjusts stepsize due to the ratios of multiple values:
%  -- needed number of iterations and optimal number of iterations
%  -- change of curvature of path and optimal change of curvature (1)
%  -- speed of continuation and optimal speed
%  -- rate of contraction and optimal rate
%  -- distance of predictor to solution point and optimal distance
%
%  The quotiens are weighted by weights specified in 'weightsMultiplicative'
%  and then multiplied.
%  The optimal values are specified by:
%  -- optimal number of iterations: 'nIterOpt'
%  -- optimal speed: 'speedOfContinuation'
%  -- rate of contraction: 'optimalContractionRate'
%  -- distance of predictor: 'predictorDistance'
%
%
%   Inputs:
%       Solver.output -- contains information of solver, such as the 
%                        needed number of iterations and rate of contraction.
%       Path          -- contains the solution points of the path and the
%                        predictors.
%       Opt           -- contains user inputs, such as values and the
%                        weights, accessible by 'weightsMultiplicative'.
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
function [xi] = multiplicativeAlt(Solver,Path,Opt)
    % create vector with Path.varAll and Path.lAll
    %
    xAll = [Path.varAll;Path.lAll];
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
    %% Weigths
    %
    weights = Opt.weightsMultiplicative;
    %
    %% Factor by number of iterations
    % correct number of iterations
    if weights(1) ~= 0 && ~isempty(Solver.output.iterations)
        if Opt.dsMax==inf
            wIter = max(Solver.output.iterations(end),1);
        else
            wIter = Solver.output.iterations(end);
        end
    else
        wIter = wTarget(1);
    end
    %
    %% Factor by contraction rate
    %
    if weights(2) ~= 0 && ~isempty(Solver.output.rateOfContraction)
        wContr = Solver.output.rateOfContraction(end);
    else
        wContr = wTarget(2);
    end
    %
    %% Factor by speed of continuation
    %
    if weights(3) ~= 0 && ~isempty(Path.speedOfContinuation)
        wSpeed = Path.speedOfContinuation(end);
    else
        wSpeed = wTarget(3);
    end
    %
    %% Factor by change of angle
    %
    % Check if there are enough solution points
    %
    if weights(4) ~= 0 && length(Path.lAll) > 3 
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
        % calculate change of angle
        %
        if angle1 <= Opt.solverTol || angle2 <= Opt.solverTol
            wCurv = 1;
        else
            wCurv = angle1/angle2;
        end
    else
        wCurv = wTarget(4);
    end
    %
    %% Factor by distance of predictor
    %
    if weights(5) ~= 0 && ~isempty(Path.xPredictor)
        if norm(xAll(:,end)) > 1e-6
            relDistanceOfPredictor = norm(Path.xPredictor(:,end) - xAll(:,end)) / norm(xAll(:,end));
            wDist = relDistanceOfPredictor;
        elseif norm(Path.xPredictor(:,end)) > 1e-6
            relDistanceOfPredictor = norm(Path.xPredictor(:,end) - xAll(:,end)) / norm(Path.xPredictor(:,end));
            wDist = relDistanceOfPredictor;
        else
            wDist = 1;
        end
    else
        wDist = wTarget(5);
    end
    %
    %% values
    %
    wI = [wIter, wContr, wSpeed, wCurv, wDist];
    %
    %% Quotients
    %
    quods = wTarget./wI;
    %
    %% limit quotients
    %
    maxQuod = 2;
    quods(quods > maxQuod) = maxQuod;
    quods(quods < 1/maxQuod) = 1/maxQuod;
    %    
    %% adjustment factor
    %
    xi = prod(quods.^weights);
end