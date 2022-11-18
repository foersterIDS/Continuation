%% path continuation - stepSize.updateTemp
%
%   Institute of Dynamics and Vibration Research
%   Leibniz University Hannover
%   28.03.2022 - Alwin Förster
%
function [Temp] = updateTemp(Path,Solver,StepsizeOptions,Temp)
    %% save old values of stepsize information
    %
    if StepsizeOptions.iterations
        Temp.length = length(Solver.output.iterations);
        Temp.flag = 0;
        if Temp.length < 3 && Temp.length > 0
            toTake = 1:Temp.length;
        elseif Temp.length == 0
            Temp.flag = 1;
        else
            toTake = Temp.length-1:Temp.length;
        end
        if isfield(Solver.output, 'iterations') && ~isempty(Solver.output.iterations) && ~Temp.flag
            Temp.iterations = Solver.output.iterations(toTake);
        else
            Temp.iterations = [];
        end
    end
    %
    if StepsizeOptions.speedOfContinuation
        Temp.length = length(Path.speedOfContinuation);
        Temp.flag = 0;
        if Temp.length < 3 && Temp.length > 0
            toTake = 1:Temp.length;
        elseif Temp.length == 0
            Temp.flag = 1;
        else
            toTake = Temp.length-1:Temp.length;
        end
        if isfield(Path, 'speedOfContinuation') && ~isempty(Path.speedOfContinuation) && ~Temp.flag
            Temp.speedOfContinuation = Path.speedOfContinuation(toTake);
        else
            Temp.speedOfContinuation = [];
        end
    end
    %
    if StepsizeOptions.predictor
        Temp.length = size(Path.xPredictor,2);
        Temp.flag = 0;
        if Temp.length < 3 && Temp.length > 0
            toTake = 1:Temp.length;
        elseif Temp.length == 0
            Temp.flag = 1;
        else
            toTake = Temp.length-1:Temp.length;
        end
        if isfield(Path, 'xPredictor') && ~isempty(Path.xPredictor) && ~Temp.flag
            Temp.predictor = Path.xPredictor(:,toTake);
        else
            Temp.predictor = [];
        end
    end
    %
    if StepsizeOptions.rateOfContraction
        Temp.length = length(Solver.output.rateOfContraction);
        Temp.flag = 0;
        if Temp.length < 3 && Temp.length > 0
            toTake = 1:Temp.length;
        elseif Temp.length == 0
            Temp.flag = 1;
        else
            toTake = Temp.length-1:Temp.length;
        end
        if isfield(Solver.output, 'rateOfContraction') && ~isempty(Solver.output.rateOfContraction) && ~Temp.flag
            Temp.rateOfContraction = Solver.output.rateOfContraction(toTake);
        else
            Temp.rateOfContraction = [];
        end
    end
    %
end