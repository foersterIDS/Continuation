%% path continuation - step_size.error
%
%   Institute of Dynamics and Vibration Research
%   Leibniz University Hannover
%   20.01.2022 - Tido Kubatschek
%
%
function [xi] = error(solver_output,Path,Opt)
    E_i = calc_error(solver_output,Path,Opt);
    
    %% adjustment factor
    E_max = 20;
    xi = 2^(E_i/E_max);
end

function E_i = calc_error(solver_output,Path,Opt)
    % create vector with Path.var_all and Path.l_all
    %
    x_all = [Path.var_all;Path.l_all];
    %
    %% define target values
    %
    % w_target: optimal number of iterations, optimal change of curvature,
    %           optimal speed of continuation, optimal rate of contraction,
    %           optimal distance of predictor
    %
    w_target = [Opt.n_iter_opt, 1, 10, 0.25, 0.5];
    %
    %% Factor by number of iterations
    % correct number of iterations
    if Opt.ds_max==inf
        w_iter = max(solver_output.iterations,1);
    else
        w_iter = solver_output.iterations;
    end    
    %% Factor by change of curvature
    %
    % Check if there are enough solution points
    %
    if length(Path.l_all) > 3
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
    w_speed = Path.speed_of_continuation(end);
    %
    %% Factor by contraction rate
    w_contr = solver_output.rate_of_contraction(end);
    %
    %% Factor by distance of predictor
    rel_distance_of_predictor = norm(Path.x_predictor(:,end) - x_all(:,end)) / norm(x_all(:,end));
    w_dist = rel_distance_of_predictor;
    %
    %% values
    w_i = [w_iter, w_curv, w_speed, w_contr, w_dist];
    %
    %% errors
    e_i = (w_target - w_i)./w_target;    
    %
    %% Weigths
    W = [0.7, 0.1, 0.1, 0.05, 0.05].';
    %
    %% weighted error
    E_i = e_i * W;
end