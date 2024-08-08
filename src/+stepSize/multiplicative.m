%% path continuation - stepSize.multiplicative
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
function [xi] = multiplicative(Solver,Path,Opt)
    % create vector with Path.varAll and Path.lAll
    %
    xAll = [Path.varAll;Path.lAll];
    %
    %% define target values
    %
    % wTarget: optimal number of iterations, optimal rate of contraction,
    %          optimal speed of continuation, optimal change of curvature,
    %          optimal distance of predictor
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
    %% Factor by change of curvature
    %
    % Check if there are enough solution points
    %
    if weights(4) ~= 0 && Path.nAll > 3 
        % calc second order derivatives of current and last step with 
        % respect to arclength
        %
        deltaS = Path.sAll(end:-1:end-2) - Path.sAll(end-1:-1:end-3);
        %
        rPp1 = 1/(deltaS(1))^2 * (xAll(:,end) - (1 + deltaS(1)/deltaS(2))*xAll(:,end-1) + deltaS(1)/deltaS(2) * xAll(:,end-2));
        rPp2 = 1/(deltaS(2))^2 * (xAll(:,end-1) - (1 + deltaS(2)/deltaS(3))*xAll(:,end-2) + deltaS(2)/deltaS(3) * xAll(:,end-3));
        %
        % calculate curvature as 2-norm
        %
        kappaCurrent = norm(rPp1);
        kappaPrevious = norm(rPp2);
        %
        % calculate change of curvature
        %
        if kappaPrevious <= Opt.solverTol || kappaCurrent <= Opt.solverTol
            wCurv = 1;
        else
            wCurv = kappaCurrent/kappaPrevious;
        end
    else
        wCurv = wTarget(4);
    end
    %
    %% Factor by distance of predictor
    %
    if weights(5) ~= 0 && ~isempty(Path.xPredictorAll)
        if norm(xAll(:,end)) > 1e-6
            relDistanceOfPredictor = norm(Path.xPredictorAll(:,end) - xAll(:,end)) / norm(xAll(:,end));
            wDist = relDistanceOfPredictor;
        elseif norm(Path.xPredictorAll(:,end)) > 1e-6
            relDistanceOfPredictor = norm(Path.xPredictorAll(:,end) - xAll(:,end)) / norm(Path.xPredictorAll(:,end));
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