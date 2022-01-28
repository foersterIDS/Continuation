%% path continuation - step_size.multiplicative
%  Adjusts stepsize due to the ratios of multiple values:
%  -- needed number of iterations and optimal number of iterations
%  -- change of curvature of path and optimal change of curvature (1)
%  -- speed of continuation and optimal speed
%  -- rate of contraction and optimal rate
%  -- distance of predictor to solution point and optimal distance
%
%  The quotiens are weighted by weights specified in 'weights_multiplicative'
%  and then multiplied.
%  The optimal values are specified by:
%  -- optimal number of iterations: 'n_iter_opt'
%  -- optimal speed: 'speed_of_continuation'
%  -- rate of contraction: 'optimal_contraction_rate'
%  -- distance of predictor: 'predictor_distance'
%
%
%   Inputs:
%       solver_output -- contains information of solver, such as the 
%                        needed number of iterations and rate of contraction.
%       Path          -- contains the solution points of the path and the
%                        predictors.
%       Opt           -- contains user inputs, such as values and the
%                        weights, accessible by 'weights_multiplicative'.
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
%   17.01.2022 - Tido Kubatschek
%
function [xi] = multiplicative(solver_output,Path,Opt)
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
    w_target = [Opt.n_iter_opt, 1, Opt.speed_of_continuation,...
        Opt.optimal_contraction_rate, Opt.predictor_distance];
    %
    %% Weigths
    %
    weights = Opt.weights_multiplicative;
    %
    %% Factor by number of iterations
    % correct number of iterations
    if weights(1) ~= 0
        if Opt.ds_max==inf
            w_iter = max(solver_output.iterations(end),1);
        else
            w_iter = solver_output.iterations(end);
        end
    else
        w_iter = w_target(1);
    end
    %
    %% Factor by change of curvature
    %
    % Check if there are enough solution points
    %
    if weights(2) ~= 0 && length(Path.l_all) > 3 
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
        kappa_current = norm(r_pp1);
        kappa_previous = norm(r_pp2);
        %
        % calculate change of curvature
        %
        if kappa_previous <= Opt.solver_tol || kappa_current <= Opt.solver_tol
            w_curv = 1;
        else
            w_curv = kappa_previous/kappa_current;
        end
    else
        w_curv = w_target(2);
    end
    %% Factor by speed of continuation
    %
    if weights(3) ~= 0
        w_speed = Path.speed_of_continuation(end);
    else
        w_speed = w_target(3);
    end
    %
    %% Factor by contraction rate
    %
    if weights(4) ~= 0
        w_contr = solver_output.rate_of_contraction(end);
    else
        w_contr = w_target(4);
    end
    %
    %% Factor by distance of predictor
    %
    if weights(5) ~= 0
        if norm(x_all(:,end)) > 1e-6
            rel_distance_of_predictor = norm(Path.x_predictor(:,end) - x_all(:,end)) / norm(x_all(:,end));
            w_dist = rel_distance_of_predictor;
        elseif norm(Path.x_predictor(:,end)) > 1e-6
            rel_distance_of_predictor = norm(Path.x_predictor(:,end) - x_all(:,end)) / norm(Path.x_predictor(:,end));
            w_dist = rel_distance_of_predictor;
        else
            w_dist = 1;
        end
    else
        w_dist = w_target(5);
    end
    %
    %% values
    %
    w_i = [w_iter, w_curv, w_speed, w_contr, w_dist];
    %
    %% Quotients
    %
    quods = w_target./w_i;
    quods(2) = w_curv;
    %
    %% limit quotients
    %
    max_quod = 2;
    quods(quods > max_quod) = max_quod;
    quods(quods < 1/max_quod) = 1/max_quod;
    %    
    %% adjustment factor
    %
    xi = prod(quods.^weights);
end