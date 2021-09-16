%% path continuation - step_size_control_curvature
%
%   Institute of Dynamics and Vibration Research
%   Leibniz University Hannover
%   21.10.2020 - Tido Kubatschek
%
function [dsn] = step_size_control_curvature(ds,ds0,Counter,solver_output,Do,Path,Opt)
    % create vector with Path.var_all and Path.l_all
    x_all = [Path.var_all;Path.l_all];
    
    % calc derivatives with respect to s
    delta_s = Path.s_all(end) - Path.s_all(end-1);
    r_p = (x_all(:,end) - x_all(:,end-1)) / delta_s;
    r_pp = (x_all(:,end-2) - 2 * x_all(:,end-1) + x_all(:,end)) / delta_s^2;
    
    % calc curvature vector
    e_1 = r_p / norm(r_p);
    e_2 = r_pp - dot(r_pp, e_1) * e_1;
    e_2 = e_2 / norm(e_2);
    
    % find max value and mean value (of abs vals) of e_2 to detect curvature
    max_curv = max(e_2);
    mean_curv = sum(abs(e_2)) / length(Path.l_all);
    
    % if value is close to zero, set it to zero
    if max_curv <= 1e-6
        max_curv = 0;
    end
    
    if mean_curv <= 1e-6
        mean_curv = 0;
    end
    
    % calc curv by weighing max and mean value
    a = 0.05;
    b = 1-a;
    curv = (a*mean_curv + b*max_curv);
    
    % calc ratio
    alpha_max = Opt.step_size_alpha_max;
    alpha = curv / alpha_max;
    
    if alpha < 1 % if ratio < 1, set alpha to 1
            alpha = 1;
    end
        
    % calc new stepsize
    dsn = ds*sqrt(Opt.n_iter_opt/(solver_output.iterations*alpha));
    dsn = max(ds/2,dsn);
    dsn = min(ds*2,dsn);    
end