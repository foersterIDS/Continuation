%% path continuation - aux.updateStepsizeData
%
%   Institute of Dynamics and Vibration Research
%   Leibniz University Hannover
%   01.02.2022 - Tido Kubatschek
%
function [Solver,Path] = updateStepsizeData(StepsizeOptions,Temp,Solver,Path)
    if StepsizeOptions.iterations
        Solver.output.iterations = Temp.iterations;
    end
    %
    if StepsizeOptions.speedOfContinuation
        Path.speedOfContinuation = Temp.speedOfContinuation;
    end
    %
    if StepsizeOptions.predictor
        Path.xPredictor = Temp.predictor;
    end
    %
    if StepsizeOptions.rateOfContraction
        Solver.output.rateOfContraction = Temp.rateOfContraction;
    end
end