%% path continuation - stepSize.updatePathStepsizeInfo
%
%   Institute of Dynamics and Vibration Research
%   Leibniz University Hannover
%   28.03.2022 - Alwin Förster
%
function [PathStepsizeInfo] = updatePathStepsizeInfo(Path,Solver,StepsizeOptions,PathStepsizeInfo)
    %% save old values of stepsize information
    %
    if StepsizeOptions.iterations
        lengthTemp = length(Solver.output.iterations);
        flagTemp = 0;
        if lengthTemp < 3 && lengthTemp > 0
            toTake = 1:lengthTemp;
        elseif lengthTemp == 0
            flagTemp = 1;
        else
            toTake = lengthTemp-1:lengthTemp;
        end
        if isfield(Solver.output, 'iterations') && ~isempty(Solver.output.iterations) && ~flagTemp
            PathStepsizeInfo.iterations = Solver.output.iterations(toTake);
        else
            PathStepsizeInfo.iterations = [];
        end
    end
    %
    if StepsizeOptions.speedOfContinuation
        lengthTemp = length(Path.speedOfContinuation);
        flagTemp = 0;
        if lengthTemp < 3 && lengthTemp > 0
            toTake = 1:lengthTemp;
        elseif lengthTemp == 0
            flagTemp = 1;
        else
            toTake = lengthTemp-1:lengthTemp;
        end
        if isfield(Path, 'speedOfContinuation') && ~isempty(Path.speedOfContinuation) && ~flagTemp
            PathStepsizeInfo.speedOfContinuation = Path.speedOfContinuation(toTake);
        else
            PathStepsizeInfo.speedOfContinuation = [];
        end
    end
    %
    if StepsizeOptions.predictor
        lengthTemp = size(Path.xPredictor,2);
        flagTemp = 0;
        if lengthTemp < 3 && lengthTemp > 0
            toTake = 1:lengthTemp;
        elseif lengthTemp == 0
            flagTemp = 1;
        else
            toTake = lengthTemp-1:lengthTemp;
        end
        if isfield(Path, 'xPredictor') && ~isempty(Path.xPredictor) && ~flagTemp
            PathStepsizeInfo.predictor = Path.xPredictor(:,toTake);
        else
            PathStepsizeInfo.predictor = [];
        end
    end
    %
    if StepsizeOptions.rateOfContraction
        lengthTemp = length(Solver.output.rateOfContraction);
        flagTemp = 0;
        if lengthTemp < 3 && lengthTemp > 0
            toTake = 1:lengthTemp;
        elseif lengthTemp == 0
            flagTemp = 1;
        else
            toTake = lengthTemp-1:lengthTemp;
        end
        if isfield(Solver.output, 'rateOfContraction') && ~isempty(Solver.output.rateOfContraction) && ~flagTemp
            PathStepsizeInfo.rateOfContraction = Solver.output.rateOfContraction(toTake);
        else
            PathStepsizeInfo.rateOfContraction = [];
        end
    end
    %
end