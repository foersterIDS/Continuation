%% path continuation - step_size_control
%
%   Institute of Dynamics and Vibration Research
%   Leibniz University Hannover
%   08.05.2020 - Alwin Förster
%   03.06.2020 - Niklas Marhenke
%   21.10.2020 - Tido Kubatschek
%
function [dsn] = step_size_control(ds,ds0,Counter,solver_output,Do,x_plus,Path,Opt)
    if ~Do.stepback
        if ~Do.deflate
            if Counter.error == 0
                if Opt.step_size_control.iterations
                    dsn = step_size_control_iterations(ds,ds0,Counter,solver_output,Do,Path,Opt);
                elseif Opt.step_size_control.angle
                    dsn = step_size_control_angle(ds,ds0,Counter,solver_output,Do,Path,Opt);
                elseif Opt.step_size_control.curvature
                    if length(Path.l_all) > 3
                        dsn = step_size_control_curvature(ds,ds0,Counter,solver_output,Do,Path,Opt);
                    else
                        dsn = step_size_control_iterations(ds,ds0,Counter,solver_output,Do,Path,Opt);
                    end
                elseif Opt.step_size_control.pid
                    if length(Path.l_all) > 4
                        dsn = step_size_control_pid(ds,ds0,Counter,solver_output,Do,Path,Opt);
                    else
                        dsn = step_size_control_iterations(ds,ds0,Counter,solver_output,Do,Path,Opt);
                    end
                elseif Opt.step_size_control.fix
                    dsn = ds0;
                else
                    error('Invalid settings for step size control!');
                end
            else
                if Opt.step_size_control.fix
                    dsn = ds0;
                else
                    dsn = ds/2;
                end
            end
        else
            dsn = ds;
        end
        %% additional penalty
        if Do.remove
            dsn = dsn*rand;
        end
        %% Limit to max./min. step size:
        dsn = min([norm(Opt.ds_max),dsn]);
        dsn = max([Opt.ds_min,dsn]);
    else
        xe = [Path.var_all(:,end);Path.l_all(end)];
        dsn = norm(x_plus-xe)/2;
    end
end