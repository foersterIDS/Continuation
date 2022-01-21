%% path continuation - step_size.angle_custom
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
    deviation_of_iterations = Opt.n_iter_opt/solver_output.iterations;
    %
    % calculate ratio
    ratio = angle / Opt.step_size_angle;
    %
    if ratio < 0.9 % if ratio < 0.9, set ch to 0.9
        ratio = 0.9;
    elseif ratio > 1.5 % high ratios are punished stronger (but ratio_max = 5)
        ratio = 5;
    end
    %
    xi = (1/ratio)^0.3 * deviation_of_iterations^0.5;
    %
end