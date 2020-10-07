%% path continuation - arclength
%
%   Institute of Dynamics and Vibration Research
%   Leibniz University Hannover
%   08.05.2020 - Alwin Förster
%   03.06.2020 - Niklas Marhenke
%
function [dsn] = arclength(ds,ds0,error_counter,solver_output,do_deflate,Opt)
    if ~do_deflate
        if error_counter==0
            dsn = ds*Opt.n_iter_opt/solver_output.iterations;
            dsn = max(ds/2,dsn);
            dsn = min(2*ds,dsn);
        else
            dsn = ds/2;
        end
    else
        dsn = ds;
    end
    %% Limit to max./min. step size:
    dsn = min([Opt.ds_max,dsn]);
    dsn = max([Opt.ds_min,dsn]);
end