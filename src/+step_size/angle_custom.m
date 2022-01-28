%% path continuation - step_size.angle_custom
%  Adjusts stepsize by the ratio of the angle of the lines connecting 
%  three consecutive solution points and an optimal angle specified in 
%  'step_size_angle'. Also adapts due to needed number of iterations 
%  (see <a href="matlab:doc('step_size.iterations_polynomial')">step_size.iterations_polynomial</a>).
%  Both adaption factors are weighted by the weights specified in 
%  'weights_angle_custom' and then multiplied.
%
%
%   Inputs:
%       solver_output -- contains information of solver, such as the 
%                        needed number of iterations.
%       Path          -- contains the solution points of the path
%       Opt           -- contains user inputs, such as optimal number of
%                        iterations, accessible by 'n_iter_opt' and the
%                        weights, accessible by 'weights_angle_custom'.
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
function [xi] = angle_custom(solver_output,Path,Opt)
    %
    % collect needed data of Path
    %
    v_needed = Path.var_all(:,end-2:end);
    l_needed = Path.l_all(end-2:end);
    z_needed = [v_needed;l_needed];
    %
    % calculate connecting vectors
    %
    v1 = z_needed(:,end) - z_needed(:,end-1);
    v2 = z_needed(:,end-1) - z_needed(:,end-2);
    %
    % calculate angle
    %
    angle = aux.vector_angle(v1,v2);
    %
    % calculate deviation of iterations
    %
    deviation_of_iterations = Opt.n_iter_opt/solver_output.iterations(end);
    %
    % calculate ratio
    ratio = angle / Opt.step_size_angle;
    %
    if ratio < 0.9 % if ratio < 0.9, set ratio to 0.9
        ratio = 0.9;
    elseif ratio > 1.2 % high ratios are punished stronger (but ratio_max = 5)
        ratio = 5;
    end
    %
    % get weigths
    weights = Opt.weights_angle_custom;
    %
    xi = deviation_of_iterations^weights(1) * (1/ratio)^weights(2);
    %
end