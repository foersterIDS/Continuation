%% path continuation - step_size.szyszkowski
%  Adjusts stepsize by the ratio of the angles of the lines connecting the
%  last four consecutive solution points. Also adapts due to needed number 
%  of iterations (see <a href="matlab:doc('step_size.iterations_polynomial')">step_size.iterations_polynomial</a>).
%  Both adaption factors are weighted by the weights specified in 
%  'weights_szyszkowski' and then multiplied.
%
%
%   Inputs:
%       solver_output -- contains information of solver, such as the 
%                        needed number of iterations.
%       Path          -- contains the solution points of the path
%       Opt           -- contains user inputs, such as optimal number of
%                        iterations, accessible by 'n_iter_opt' and the
%                        weights, accessible by 'weights_szyszkowski'.
%                        
%   Outputs:
%       xi            -- stepsize adaption factor
%
%
%
%  See the <a href="matlab:open('..\doc\html\continuation.html')">documentation</a>. See <a href="matlab:doc('step_size.control')">other stepsize adaption methods</a>.
%
%   DOI: https://doi.org/10.1007/s004660050513  (adapted version)
%
%   Institute of Dynamics and Vibration Research
%   Leibniz University Hannover
%   17.01.2022 - Tido Kubatschek
%
function [xi] = szyszkowski(solver_output,Path,Opt)
    %% Method of Szyszkowski and Husband
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
    angle = aux.vector_angle(v1,v2);
    angle_m1 = aux.vector_angle(v2,v3);
    %
    % calculate next angle
    %
    angle_p1 = 2*angle - angle_m1 * norm(v1)/norm(v2);
    %
    % check if angle is too large
    %
    if angle > Opt.step_size_angle || angle_p1 > Opt.step_size_angle
        xi = 0.5;
    else
        % correct number of iterations
        %
        if Opt.ds_max==inf
            iter = max(solver_output.iterations(end),1);
        else
            iter = solver_output.iterations(end);
        end
        %
        % calculate deviation of iterations
        %
        deviation_of_iterations = Opt.n_iter_opt/iter;
        %
        % get weigths
        weights = Opt.weights_szyszkowski;
        %
        % calculate new step size
        %
        xi = deviation_of_iterations^weights(1) *...
            abs(angle/angle_p1)^weights(2);
        %
    end
end