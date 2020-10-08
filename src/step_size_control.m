%% path continuation - step_size_control
%
%   Institute of Dynamics and Vibration Research
%   Leibniz University Hannover
%   08.05.2020 - Alwin Förster
%   03.06.2020 - Niklas Marhenke
%
function [dsn] = step_size_control(ds,ds0,error_counter,solver_output,do_deflate,vars,ls,Opt)
    if ~do_deflate
        if error_counter == 0
            if Opt.step_size_control.standard
                dsn = step_size_control_standard(ds,ds0,error_counter,solver_output,do_deflate,vars,ls,Opt);
            elseif Opt.step_size_control.angle
                dsn = step_size_control_angle(ds,ds0,error_counter,solver_output,do_deflate,vars,ls,Opt);
            else
                error('Invalid settings for step size control!');
            end
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