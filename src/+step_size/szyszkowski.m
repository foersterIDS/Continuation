%% path continuation - step_size.szyszkowski
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
    % calculate next angle
    %
    angle_p1 = 2*angle - angle_m1 * norm(v1)/norm(v2);
    %
    % check if angle is too large
    %
    if angle > Opt.step_size_angle
        xi = 1/2;
    else
        % calculate deviation of iterations
        %
        deviation_of_iterations = Opt.n_iter_opt/solver_output.iterations;
        %
        % calculate new step size
        %
        xi = abs(angle/angle_p1) * deviation_of_iterations^0.5;
        %
    end
end