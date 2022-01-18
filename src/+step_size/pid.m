%% path continuation - step_size.pid
%
%   Institute of Dynamics and Vibration Research
%   Leibniz University Hannover
%   21.10.2020 - Tido Kubatschek
%
%   DOI: 10.1002/ﬂd.1998
%
function [dsn] = pid(ds,ds0,Counter,solver_output,Do,Path,Opt)
    % parameters of pid control 
    k_P = 0.1;
    k_I = 0.01;
    k_D = 0.5;
    
    pid_tol = Opt.step_size_pid_tol;
    
    dvarsdl = @(k) (Path.var_all(:,end+k) - Path.var_all(:,end+k-1)) / (Path.l_all(end+k) - Path.l_all(end+k-1));
    
    e_star = @(k) norm(dvarsdl(k) - dvarsdl(k-1)) / norm(dvarsdl(k));
    
    e = [e_star(0), e_star(-1), e_star(-2)];

    e = e / pid_tol;
    
    % calculate step size
    dsn = (e(2) / e(1))^k_P * (1 / e(1))^k_I * (e(2)^2 / (e(1) * e(3)) )^k_D * ds;
end