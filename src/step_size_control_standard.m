%% path continuation - step_size_control_standard
%
%   Institute of Dynamics and Vibration Research
%   Leibniz University Hannover
%   08.10.2020 - Alwin Förster
%
function [dsn] = step_size_control_standard(ds,ds0,error_counter,solver_output,do_deflate,vars,ls,Opt)
    % calculate step size
    dsn = ds*Opt.n_iter_opt/(solver_output.iterations);
    dsn = max(ds/2,dsn);
    dsn = min(ds*2,dsn);
end