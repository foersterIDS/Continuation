%% path continuation - step_size.angle_change
%  Adjusts stepsize by the change of the angles of the lines connecting the
%  last four consecutive solution points. Also adapts due to needed number 
%  of iterations (see <a href="matlab:doc('step_size.iterations_polynomial')">step_size.iterations_polynomial</a>).
%  Both adaption factors are weighted by the weights specified in 
%  'weights_angle_change' and then multiplied.
%
%
%   Inputs:
%       Solver.output -- contains information of solver, such as the 
%                        needed number of iterations.
%       Path          -- contains the solution points of the path
%       Opt           -- contains user inputs, such as optimal number of
%                        iterations, accessible by 'n_iter_opt' and the
%                        weights, accessible by 'weights_angle_change'.
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
function [xi] = angle_change(Solver,Path,Opt)
    %
    % collect needed data of Path
    %
    var_needed = Path.var_all(:,end-3:end);
    l_needed = Path.l_all(end-3:end);
    z_needed = [var_needed; l_needed];
    %
    % calculate connecting vectors
    %
    v1 = z_needed(:,end) - z_needed(:,end-1);
    v2 = z_needed(:,end-1) - z_needed(:,end-2);
    v3 = z_needed(:,end-2) - z_needed(:,end-3);
    %
    % calculate angles
    %
    angle_1 = aux.vector_angle(v1,v2);
    angle_2 = aux.vector_angle(v2,v3);
    %
    % calculate change of curvature
    %
    change_of_angle = angle_2/angle_1;
    %
    % correct number of iterations
    %
    if Opt.ds_max==inf
        iter = max(Solver.output.iterations(end),1);
    else
        iter = Solver.output.iterations(end);
    end
    %
    % calculate deviation of iterations
    %
    deviation_of_iterations = Opt.n_iter_opt/iter;
    %
    % get weigths
    %
    weights = Opt.weights_angle_change;
    %
    % calculate adaption factor
    %
    xi = deviation_of_iterations^weights(1) * change_of_angle^weights(2);
    %
end