%% path continuation - step_size.fayezioghani
%
%   Institute of Dynamics and Vibration Research
%   Leibniz University Hannover
%   17.01.2022 - Tido Kubatschek
%
function [xi] = fayezioghani(ds,solver_output,Path,Jac,Opt)
    %% Method of Fayezioghani et al.
    %
    % collect needed data of Path
    %
    v_needed = Path.var_all(:,end-2:end);
    l_needed = Path.l_all(end-2:end);
    z_needed = [v_needed; l_needed];
    %
    % calculate connecting vector and tangent
    %
    v = z_needed(:,end) - z_needed(:,end-1);
    [~,tangent] = predictor_ode(Path,ds,Jac,[]);
    %
    % calculate angle
    %
    angle = vector_angle(v,tangent);
    %
    % calculate new step size
    %
    beta_1 = 0.5;
    beta_2 = 0.1;
    xi = (Opt.n_iter_opt/(solver_output.iterations))^beta_1 *...
        ((cos(angle) +1)/(cos(Opt.step_size_angle) +1))^beta_2;
end