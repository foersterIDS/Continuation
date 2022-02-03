%% path continuation - step_size.multiplicative_alt
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
%       Solver.output -- contains information of solver, such as the 
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
function [xi] = multiplicative_alt(Solver,Path,Opt)
    % create vector with Path.var_all and Path.l_all
    %
    x_all = [Path.var_all;Path.l_all];
    %
    %% define target values
    %
    % w_target: optimal number of iterations, optimal rate of contraction,
    %           optimal speed of continuation, optimal change of angle,
    %           optimal distance of predictor
    %
    w_target = [Opt.n_iter_opt, Opt.optimal_contraction_rate,...
        Opt.speed_of_continuation, 1, Opt.predictor_distance];
    %
    %% Weigths
    %
    weights = Opt.weights_multiplicative;
    %
    %% Factor by number of iterations
    % correct number of iterations
    if weights(1) ~= 0 && ~isempty(Solver.output.iterations)
        if Opt.ds_max==inf
            w_iter = max(Solver.output.iterations(end),1);
        else
            w_iter = Solver.output.iterations(end);
        end
    else
        w_iter = w_target(1);
    end
    %
    %% Factor by contraction rate
    %
    if weights(2) ~= 0 && ~isempty(Solver.output.rate_of_contraction)
        w_contr = Solver.output.rate_of_contraction(end);
    else
        w_contr = w_target(2);
    end
    %
    %% Factor by speed of continuation
    %
    if weights(3) ~= 0 && ~isempty(Path.speed_of_continuation)
        w_speed = Path.speed_of_continuation(end);
    else
        w_speed = w_target(3);
    end
    %
    %% Factor by change of angle
    %
    % Check if there are enough solution points
    %
    if weights(4) ~= 0 && length(Path.l_all) > 3 
        %
        % calculate connecting vectors
        %
        v1 = x_all(:,end) - x_all(:,end-1);
        v2 = x_all(:,end-1) - x_all(:,end-2);
        v3 = x_all(:,end-2) - x_all(:,end-3);
        %
        % calculate angles
        %
        angle_1 = aux.vector_angle(v1,v2);
        angle_2 = aux.vector_angle(v2,v3);
        %
        % calculate change of angle
        %
        if angle_1 <= Opt.solver_tol || angle_2 <= Opt.solver_tol
            w_curv = 1;
        else
            w_curv = angle_1/angle_2;
        end
    else
        w_curv = w_target(4);
    end
    %
    %% Factor by distance of predictor
    %
    if weights(5) ~= 0 && ~isempty(Path.x_predictor)
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
    w_i = [w_iter, w_contr, w_speed, w_curv, w_dist];
    %
    %% Quotients
    %
    quods = w_target./w_i;
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