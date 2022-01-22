%% path continuation - step_size.iterations_exponential
%
%   Institute of Dynamics and Vibration Research
%   Leibniz University Hannover
%   08.10.2020 - Alwin Förster
%   17.01.2022 - Tido Kubatschek
%
%
function [xi] = iterations_exponential(solver_output,Opt)
    % correct number of iterations
    if Opt.ds_max==inf
        iter = max(solver_output.iterations,1);
    else
        iter = solver_output.iterations;
    end
    % calculate step size adaption factor
    xi = 2^((Opt.n_iter_opt - iter)/Opt.step_size_exponential_weigth);
end