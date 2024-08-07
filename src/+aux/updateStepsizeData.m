%% path continuation - aux.updateStepsizeData
%
%   Institute of Dynamics and Vibration Research
%   Leibniz University Hannover
%   01.02.2022 - Tido Kubatschek
%
function [Solver,Path] = updateStepsizeData(StepsizeOptions,PathStepsizeInfo,Solver,Path)
    if StepsizeOptions.iterations
        Solver.output.iterations = PathStepsizeInfo.iterations;
    end
    %
    if StepsizeOptions.speedOfContinuation
        Path.speedOfContinuation = PathStepsizeInfo.speedOfContinuation;
    end
    %
    if StepsizeOptions.predictor
        Path.xPredictor = PathStepsizeInfo.predictor;
    end
    %
    if StepsizeOptions.rateOfContraction
        Solver.output.rateOfContraction = PathStepsizeInfo.rateOfContraction;
    end
end