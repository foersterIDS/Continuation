%% path continuation - step_size_control_iterations
%
%   Institute of Dynamics and Vibration Research
%   Leibniz University Hannover
%   03.06.2020 - Niklas Marhenke
%   08.10.2020 - Alwin Förster
%
%   DOI: 10.1007/978-1-4419-1740-9
%
function [dsn] = step_size_control_iterations(ds,ds0,error_counter,solver_output,do_deflate,vars,ls,Opt)
    % calculate step size
    dsn = ds*Opt.n_iter_opt/(solver_output.iterations);
    dsn = max(ds/2,dsn);
    dsn = min(ds*2,dsn);
end