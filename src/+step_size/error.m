%% path continuation - step_size.error
%
%   Institute of Dynamics and Vibration Research
%   Leibniz University Hannover
%   20.01.2022 - Tido Kubatschek
%
%
function [xi] = error(solver_output,Path,Opt)
    %
    E_i = calc_error(solver_output,Path,Opt,0);
    %
    %
    K = Opt.step_size_error_pd;
    %
    if length(Path.l_all) > 3 && K(2) > 0
        E_i_m1 = calc_error(solver_output,Path,Opt,1);
    else
        E_i_m1 = 0;
    end
    %
    E = K(1) * E_i + K(2) * (E_i - E_i_m1);
    %
    %% adjustment factor
    %
    xi = 2^(E/Opt.step_size_error_max);
end

function E_i = calc_error(solver_output,Path,Opt,previous)
    %
    % determine length
    length_arrays = length(Path.l_all);
    %
    % determine state
    if previous
        end_of_array = 1;
        length_arrays = length_arrays - 1;
    else
        if length_arrays < 3
            end_of_array = 1;
        else
            end_of_array = 2;
        end
    end
    
    % create vector with Path.var_all and Path.l_all
    %
    x_all = [Path.var_all(:,1:length_arrays);Path.l_all(1:length_arrays)];
    %
    %% define target values
    %
    % w_target: optimal number of iterations, optimal change of curvature,
    %           optimal speed of continuation, optimal rate of contraction,
    %           optimal distance of predictor
    %
    w_target = [Opt.n_iter_opt, 1, Opt.speed_of_continuation,...
        Opt.optimal_contraction_rate, Opt.predictor_distance];
    %
    %% Weights
    weights = Opt.weights_error.';
    %
    %% Factor by number of iterations
    % correct number of iterations
    if weights(1) ~= 0
        if Opt.ds_max==inf
            w_iter = max(solver_output.iterations(end_of_array),1);
        else
            w_iter = solver_output.iterations(end_of_array);
        end
    else
        w_iter = w_target(1);
    end
    %% Factor by change of curvature
    %
    % Check if there are enough solution points
    %
    if weights(2) ~= 0 && length_arrays > 3
        % calc second order derivatives of current and last step with 
        % respect to arclength
        %
        delta_s = Path.s_all(end:-1:end-2) - Path.s_all(end-1:-1:end-3);
        %
        r_pp1 = 1/(delta_s(1))^2 * (x_all(:,end) - (1 + delta_s(1)/delta_s(2))*x_all(:,end-1) + delta_s(1)/delta_s(2) * x_all(:,end-2));
        r_pp2 = 1/(delta_s(2))^2 * (x_all(:,end-1) - (1 + delta_s(2)/delta_s(3))*x_all(:,end-2) + delta_s(2)/delta_s(3) * x_all(:,end-3));
        %
        % calculate curvature as 2-norm
        %
        kappa_current = norm(r_pp1) / norm(x_all(:,end));
        kappa_previous = norm(r_pp2) / norm(x_all(:,end-1));
        %
        % calculate change of curvature
        %
        if kappa_previous <= Opt.solver_tol || kappa_current <= Opt.solver_tol
            w_curv = 1;
        else
            w_curv = kappa_previous/kappa_current;
        end
    else
        w_curv = 1;
    end    
    %% Factor by speed of continuation
    if weights(3) ~= 0
        w_speed = Path.speed_of_continuation(end_of_array);
    else
        w_speed = w_target(3);
    end
    %
    %% Factor by contraction rate
    if weights(4) ~= 0
        w_contr = solver_output.rate_of_contraction(end_of_array);
    else
        w_contr = w_target(4);
    end
    %
    %% Factor by distance of predictor
    if weights(5) ~= 0
        if norm(x_all(:,end_of_array)) > 1e-6
            rel_distance_of_predictor = norm(Path.x_predictor(:,end_of_array) - x_all(:,end_of_array)) / norm(x_all(:,end_of_array));
            w_dist = rel_distance_of_predictor;
        elseif norm(Path.x_predictor(:,end_of_array)) > 1e-6
            rel_distance_of_predictor = norm(Path.x_predictor(:,end_of_array) - x_all(:,end_of_array)) / norm(Path.x_predictor(:,end_of_array));
            w_dist = rel_distance_of_predictor;
        else
            w_dist = 1;
        end
    else
        w_dist = w_target(5);
    end
    %
    %% values
    w_i = [w_iter, w_curv, w_speed, w_contr, w_dist];
    %
    %% errors
    e_i = (w_target - w_i)./w_target;    
    %
    %% weighted error
    E_i = (e_i * weights)/sum(weights);
end