%% path continuation - step_size.yoon
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
    angle = vector_angle(v1,v2);
    angle_m1 = vector_angle(v2,v3);
    %
    % adapt stepsize
    %
    xi_u = 1.1; xi_l = 1e10;
    a = (xi_u - 1) / (1 - xi_l);
    xi = (xi_u - ((xi_u - xi_l) * a) / ((angle_m1/angle) + a));
    %
    % calculate deviation of iterations
    %
    deviation_of_iterations = Opt.n_iter_opt/solver_output.iterations;
    %
    % correct stepsize by iterations
    %
    xi = xi * deviation_of_iterations^0.5;
    %
end