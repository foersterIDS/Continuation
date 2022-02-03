%% path continuation - step_size.error
%  Adjusts stepsize due to the relative differences of multiple values:
%  -- needed number of iterations and optimal number of iterations
%  -- change of curvature of path and optimal change of curvature (1)
%  -- speed of continuation and optimal speed
%  -- rate of contraction and optimal rate
%  -- distance of predictor to solution point and optimal distance
%  The optimal values are specified by:
%  -- optimal number of iterations: 'n_iter_opt'
%  -- optimal speed: 'speed_of_continuation'
%  -- rate of contraction: 'optimal_contraction_rate'
%  -- distance of predictor: 'predictor_distance'
%  The differences are then devided by the optimal values, weighted by 
%  weights specified in 'weights_multiplicative' and then added togehter.
%  The sum is then used as an exponent and can be weighted by
%  'step_size_error_max'. Also it is possible to consider the tendency of
%  the adaption in the last step by a PD controller which constants are 
% specified in 'step_size_error_pd'.
%
%
%   Inputs:
%       Solver.output -- contains information of solver, such as the 
%                        needed number of iterations and rate of contraction.
%       Path          -- contains the solution points of the path and the
%                        predictors.
%       Opt           -- contains user inputs, such as optimal values, the
%                        weights, accessible by 'weights_multiplicative',
%                        the sum weight ('step_size_error_max') and the 
%                        PD constants ('step_size_error_pd').
%                        
%   Outputs:
%       xi            -- stepsize adaption factor
%
%
%
%  See the <a href="matlab:open('..\doc\html\continuation.html')">documentation</a>. See <a href="matlab:doc('step_size.control')">other stepsize adaption methods</a>.
%
%   Institute of Dynamics and Vibration Research
%   Leibniz University Hannover
%   20.01.2022 - Tido Kubatschek
%
function [xi] = error(Solver,Path,Opt)
    %
    E_i = calc_error(Solver,Path,Opt,0);
    %
    %
    K = Opt.step_size_error_pd;
    %
    if length(Path.l_all) > 1 && ~isempty(Path.speed_of_continuation) && K(2) > 0
        E_i_m1 = calc_error(Solver,Path,Opt,1);
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

function E_i = calc_error(Solver,Path,Opt,previous)
    %
    %% Weights
    weights = Opt.weights_error.';
    %
    %
    % determine length
    length_path = length(Path.l_all);
    length_iterations = length(Solver.output.iterations);
    if weights(3) ~= 0
        length_arrays = length(Path.speed_of_continuation);
    elseif weights(4) ~= 0
        length_arrays = length(Solver.output.rate_of_convergence);
    elseif weights(5) ~= 0
        length_arrays = length(Path.x_predictor(1,:));
    else
        length_arrays = 1;
    end
    
    %
    % determine state
    if previous
        end_of_iterations = length_iterations - 1;
        end_of_array = length_arrays - 1;
        length_path = length_path - 1;
    else
        if length_path < 3
            end_of_iterations = length_iterations - 1;
            end_of_array = length_arrays - 1;
        else
            end_of_iterations = length_iterations;
            end_of_array = length_arrays;
        end
    end
    
    % create vector with Path.var_all and Path.l_all
    %
    x_all = [Path.var_all(:,1:length_path);Path.l_all(1:length_path)];
    %
    %% define target values
    %
    % w_target: optimal number of iterations, optimal rate of contraction,
    %           optimal speed of continuation, optimal change of curvature,
    %           optimal distance of predictor
    %
    w_target = [Opt.n_iter_opt, Opt.optimal_contraction_rate,...
        Opt.speed_of_continuation, 1, Opt.predictor_distance];
    %
    %
    %% Factor by number of iterations
    % correct number of iterations
    if weights(1) ~= 0 && ~isempty(Solver.output.iterations)
        if Opt.ds_max==inf
            w_iter = max(Solver.output.iterations(end_of_iterations),1);
        else
            w_iter = Solver.output.iterations(end_of_iterations);
        end
    else
        w_iter = w_target(1);
    end
    %
    %% Factor by contraction rate
    if weights(2) ~= 0 && ~isempty(Solver.output.rate_of_contraction)
        w_contr = Solver.output.rate_of_contraction(end_of_array);
    else
        w_contr = w_target(2);
    end
    %
    %% Factor by speed of continuation
    if weights(3) ~= 0 && ~isempty(Path.speed_of_continuation)
        w_speed = Path.speed_of_continuation(end_of_array);
    else
        w_speed = w_target(3);
    end
    %
    %% Factor by change of curvature
    %
    % Check if there are enough solution points
    %
    if weights(4) ~= 0 && length_path > 3
        % calc second order derivatives of current and last step with 
        % respect to arclength
        %
        delta_s = Path.s_all(length_path:-1:length_path-2) - Path.s_all(length_path-1:-1:length_path-3);
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
            w_curv = kappa_current/kappa_previous;
        end
    else
        w_curv = 1;
    end
    %
    %% Factor by distance of predictor
    %
    if weights(5) ~= 0 && ~isempty(Path.x_predictor)
        if norm(x_all(:,length_path)) > 1e-6
            rel_distance_of_predictor = norm(Path.x_predictor(:,end_of_array) - x_all(:,length_path)) / norm(x_all(:,length_path));
            w_dist = rel_distance_of_predictor;
        elseif norm(Path.x_predictor(:,end_of_array)) > 1e-6
            rel_distance_of_predictor = norm(Path.x_predictor(:,end_of_array) - x_all(:,length_path)) / norm(Path.x_predictor(:,end_of_array));
            w_dist = rel_distance_of_predictor;
        else
            w_dist = w_target(5);
        end
    else
        w_dist = w_target(5);
    end
    %
    %% values
    %
    w_i = [w_iter, w_contr, w_speed, w_curv, w_dist];
    %
    %% errors
    %
    e_i = (w_target - w_i)./w_target;    
    %
    %% weighted error sum
    %
    E_i = (e_i * weights)/sum(weights);
end