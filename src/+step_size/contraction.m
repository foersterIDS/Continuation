%% path continuation - step_size.contraction
%
%   Institute of Dynamics and Vibration Research
%   Leibniz University Hannover
%   20.01.2022 - Tido Kubatschek
%
%
function [xi] = contraction(solver_output,Opt)
    % get rates of contraction
    %
    % current rate
    %
    current_rate = solver_output.rate_of_contraction(end);
    %
    % if there is more than one rate accessible also get previous rate
    % and calculate relative difference
    %
    if length(solver_output.rate_of_contraction) > 1
        previous_rate = solver_output.rate_of_contraction(end-1);
        rel_difference = abs(current_rate - previous_rate) / abs(previous_rate);
    else
        rel_difference = inf;
    end
    % 
    % optimal rate
    %
    q_opt = Opt.optimal_contraction_rate;
    %
    if current_rate > q_opt && rel_difference > 0.01
        % only decrease if there was a difference
        xi = 0.5;
    elseif current_rate <= q_opt/4 || (current_rate <= q_opt/2 && rel_difference < 0.01)
        % increase if rate is small or there is just a small difference
        xi = 1.2;
    else
        xi = 1;
    end
end