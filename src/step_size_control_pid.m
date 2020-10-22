%% path continuation - step_size_control_pid
%
%   Institute of Dynamics and Vibration Research
%   Leibniz University Hannover
%   21.10.2020 - Tido Kubatschek
%
function [dsn] = step_size_control_pid(ds,ds0,error_counter,solver_output,do_deflate,vars,ls,Opt)
    % parameters of pid control
    k_P = 0.075;
    k_I = 0.175;
    k_D = 0.01;
    
    e_max = Opt.step_size_e_max;
    
    dvarsdl = @(k) (vars(:,end+k) - vars(:,end+k-1)) / (ls(end+k) - ls(end+k-1));
    
    e_star = @(k) norm(dvarsdl(k) - dvarsdl(k-1)) / norm(dvarsdl(k));
    
    e = [e_star(0), e_star(-1), e_star(-2)];

    e = e / e_max;
    
    % calculate step size
    dsn = (e(2) / e(1))^k_P * (1 / e(1))^k_I * (e(2)^2 / (e(1) * e(3)) )^k_D * ds;
    dsn = max(ds/2,dsn);
    dsn = min(ds*2,dsn);
end