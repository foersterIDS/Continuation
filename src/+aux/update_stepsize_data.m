%% path continuation - aux.update_stepsize_data
%
%   Institute of Dynamics and Vibration Research
%   Leibniz University Hannover
%   01.02.2022 - Tido Kubatschek
%
function [solver_output,Path] = update_stepsize_data(Stepsize_options,Tmp_struct,solver_output,Path)
    if Stepsize_options.iterations
        solver_output.iterations = Tmp_struct.iterations_tmp;
    end
    %
    if Stepsize_options.speed_of_continuation
        Path.speed_of_continuation = Tmp_struct.speed_of_continuation_tmp;
    end
    %
    if Stepsize_options.predictor
        Path.x_predictor = Tmp_struct.predictor_tmp;
    end
    %
    if Stepsize_options.rate_of_contraction
        solver_output.rate_of_contraction = Tmp_struct.rate_of_contraction_tmp;
    end
end