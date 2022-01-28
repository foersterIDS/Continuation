%% path continuation - step_size.yoon
%  Adjusts stepsize by the ratio of the angles of the lines connecting the
%  last four consecutive solution points. Also adapts due to needed number 
%  of iterations (see <a href="matlab:doc('step_size.iterations_polynomial')">step_size.iterations_polynomial</a>).
%  Both adaption factors are weighted by the weights specified in 
%  'weights_yoon' and then multiplied.
%
%
%   Inputs:
%       solver_output -- contains information of solver, such as the 
%                        needed number of iterations.
%       Path          -- contains the solution points of the path
%       Opt           -- contains user inputs, such as optimal number of
%                        iterations, accessible by 'n_iter_opt' and the
%                        weights, accessible by 'weights_yoon'.
%                        
%   Outputs:
%       xi            -- stepsize adaption factor
%
%
%
%  See the <a href="matlab:open('..\doc\html\continuation.html')">documentation</a>. See <a href="matlab:doc('step_size.control')">other stepsize adaption methods</a>.
%
%   DOI: 10.1177/0954406215586588  (adapted version)
%
%   Institute of Dynamics and Vibration Research
%   Leibniz University Hannover
%   17.01.2022 - Tido Kubatschek
%
function [xi] = yoon(solver_output,Path,Opt)
    %% Method of Yoon and Kim
    %
    % collect needed data of Path
    %
    v_needed = Path.var_all(:,end-3:end);
    l_needed = Path.l_all(end-3:end);
    z_needed = [v_needed; l_needed];
    %
    % calculate connecting vectors
    %
    v1 = z_needed(:,end) - z_needed(:,end-1);
    v2 = z_needed(:,end-1) - z_needed(:,end-2);
    v3 = z_needed(:,end-2) - z_needed(:,end-3);
    %
    % calculate angles
    %
    angle = aux.vector_angle(v1,v2);
    angle_m1 = aux.vector_angle(v2,v3);
    %
    % adapt stepsize
    %
    xi_u = 1.1; xi_l = 1e10;
    a = (xi_u - 1) / (1 - xi_l);
    xi = (xi_u - ((xi_u - xi_l) * a) / ((angle_m1/angle) + a));
    %
    % calculate deviation of iterations
    %
    deviation_of_iterations = Opt.n_iter_opt/solver_output.iterations(end);
    %
    % correct stepsize by iterations
    %
    xi = xi * deviation_of_iterations^Opt.weights_yoon;
    %
end