%% path continuation - aux.update_stepsize_data
%
%   Institute of Dynamics and Vibration Research
%   Leibniz University Hannover
%   01.02.2022 - Tido Kubatschek
%
function [Solver,Path] = update_stepsize_data(Stepsize_options,Temp,Solver,Path)
    if Stepsize_options.iterations
        Solver.output.iterations = Temp.iterations;
    end
    %
    if Stepsize_options.speed_of_continuation
        Path.speed_of_continuation = Temp.speed_of_continuation;
    end
    %
    if Stepsize_options.predictor
        Path.x_predictor = Temp.predictor;
    end
    %
    if Stepsize_options.rate_of_contraction
        Solver.output.rate_of_contraction = Temp.rate_of_contraction;
    end
end