%% path continuation - step_size.update_temp
%
%   Institute of Dynamics and Vibration Research
%   Leibniz University Hannover
%   28.03.2022 - Alwin Förster
%
function [Temp] = update_temp(Path,Solver,Stepsize_options,Temp)
    %% save old values of stepsize information
    %
    if Stepsize_options.iterations
        Temp.length = length(Solver.output.iterations);
        Temp.flag = 0;
        if Temp.length < 3 && Temp.length > 0
            to_take = 1:Temp.length;
        elseif Temp.length == 0
            Temp.flag = 1;
        else
            to_take = Temp.length-1:Temp.length;
        end
        if isfield(Solver.output, 'iterations') && ~isempty(Solver.output.iterations) && ~Temp.flag
            Temp.iterations = Solver.output.iterations(to_take);
        else
            Temp.iterations = [];
        end
    end
    %
    if Stepsize_options.speed_of_continuation
        Temp.length = length(Path.speed_of_continuation);
        Temp.flag = 0;
        if Temp.length < 3 && Temp.length > 0
            to_take = 1:Temp.length;
        elseif Temp.length == 0
            Temp.flag = 1;
        else
            to_take = Temp.length-1:Temp.length;
        end
        if isfield(Path, 'speed_of_continuation') && ~isempty(Path.speed_of_continuation) && ~Temp.flag
            Temp.speed_of_continuation = Path.speed_of_continuation(to_take);
        else
            Temp.speed_of_continuation = [];
        end
    end
    %
    if Stepsize_options.predictor
        Temp.length = size(Path.x_predictor,2);
        Temp.flag = 0;
        if Temp.length < 3 && Temp.length > 0
            to_take = 1:Temp.length;
        elseif Temp.length == 0
            Temp.flag = 1;
        else
            to_take = Temp.length-1:Temp.length;
        end
        if isfield(Path, 'x_predictor') && ~isempty(Path.x_predictor) && ~Temp.flag
            Temp.predictor = Path.x_predictor(:,to_take);
        else
            Temp.predictor = [];
        end
    end
    %
    if Stepsize_options.rate_of_contraction
        Temp.length = length(Solver.output.rate_of_contraction);
        Temp.flag = 0;
        if Temp.length < 3 && Temp.length > 0
            to_take = 1:Temp.length;
        elseif Temp.length == 0
            Temp.flag = 1;
        else
            to_take = Temp.length-1:Temp.length;
        end
        if isfield(Solver.output, 'rate_of_contraction') && ~isempty(Solver.output.rate_of_contraction) && ~Temp.flag
            Temp.rate_of_contraction = Solver.output.rate_of_contraction(to_take);
        else
            Temp.rate_of_contraction = [];
        end
    end
    %
end