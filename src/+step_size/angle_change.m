%% path continuation - step_size.angle_change
%
%   Institute of Dynamics and Vibration Research
%   Leibniz University Hannover
%   17.01.2022 - Tido Kubatschek
%
function [xi] = angle_change(solver_output,Path,Opt)
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
    angle_1 = vector_angle(v1,v2);
    angle_2 = vector_angle(v2,v3);
    %
    % calculate change of curvature
    %
    change_of_angle = angle_2/angle_1;
    %
    % calculate deviation of iterations
    %
    deviation_of_iterations = Opt.n_iter_opt/solver_output.iterations;
    %
    % calculate adaption factor
    %
    xi = change_of_angle^0.5 * deviation_of_iterations^0.5;
    %
end