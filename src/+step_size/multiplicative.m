%% path continuation - step_size.multiplicative
%
%   Institute of Dynamics and Vibration Research
%   Leibniz University Hannover
%   17.01.2022 - Tido Kubatschek
%
function [xi] = multiplicative(solver_output,Path,Opt)
    % create vector with Path.var_all and Path.l_all
    %
    x_all = [Path.var_all;Path.l_all];
    %
    %% Factor by number of iterations
    % correct number of iterations
    if Opt.ds_max==inf
        iter = max(solver_output.iterations,1);
    else
        iter = solver_output.iterations;
    end
    % calculate step size adaption factor
    xi_iter = (Opt.n_iter_opt/iter)^Opt.step_size_iterations_beta;
    
    %% Factor by change of curvature
    %
    % Check if there are enough solution points
    %
    if length(Path.l_all) > 3
        % calc second order derivatives of current and last step with 
        % respect to arclength
        delta_s = Path.s_all(end:-1:end-2) - Path.s_all(end-1:-1:end-3);

        r_pp1 = 1/(delta_s(1))^2 * (x_all(:,end) - (1 + delta_s(1)/delta_s(2))*x_all(:,end-1) + delta_s(1)/delta_s(2) * x_all(:,end-2));
        r_pp2 = 1/(delta_s(2))^2 * (x_all(:,end-1) - (1 + delta_s(2)/delta_s(3))*x_all(:,end-2) + delta_s(2)/delta_s(3) * x_all(:,end-3));

        % calculate curvature as 2-norm
        kappa_current = norm(r_pp1);
        kappa_previous = norm(r_pp2);

        % calculate change of curvature
        xi_curv = kappa_previous/kappa_current;
    else
        xi_curv = 1;
    end    
    %% Factor by speed of continuation
    xi_speed = 0.5/Path.speed_of_continuation(end);
    
    %% Factor by contraction rate
    xi_contr = 0.25 / solver_output.rate_of_contraction;
    
    %% Factor by distance of predictor
    rel_distance_of_predictor = norm(Path.x_predictor(:,end) - x_all(:,end)) / norm(x_all(:,end));
    xi_dist = 1e-4/rel_distance_of_predictor;
   
    %% Weigths
    W = [0.5, 0.2, 0.1, 0.2, 0.1];
    
    %% adjustment factor
    xi = prod([xi_iter, xi_curv, xi_speed, xi_contr, xi_dist].^W);
end