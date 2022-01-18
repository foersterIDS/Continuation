%% path continuation - step_size.pid_valli
%
%   Institute of Dynamics and Vibration Research
%   Leibniz University Hannover
%   21.10.2020 - Tido Kubatschek
%
%   DOI: 10.1002/ﬂd.1998
%
function [xi] = pid_valli(Path,Opt)
    dvarsdl = @(k) (Path.var_all(:,end+k) - Path.var_all(:,end+k-1)) / (Path.l_all(end+k) - Path.l_all(end+k-1));
    e_star = @(k) norm(dvarsdl(k) - dvarsdl(k-1)) / norm(dvarsdl(k));
    %
    % parameters of pid control 
    %
    PID = Opt.step_size_pid_params;
    pid_tol = Opt.step_size_pid_tol;
    %
    % relative change
    %
    e = [e_star(0), e_star(-1), e_star(-2)];
    e = e / pid_tol;
    %
    % calculate step size
    %
    xi = (e(2) / e(1))^PID(1) * (1 / e(1))^PID(2) * (e(2)^2 / (e(1) * e(3)) )^PID(3);
    %
end