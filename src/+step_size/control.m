%% path continuation - step_size.control
%
%   Institute of Dynamics and Vibration Research
%   Leibniz University Hannover
%   08.05.2020 - Alwin Förster
%   03.06.2020 - Niklas Marhenke
%   21.10.2020 - Tido Kubatschek
%
function [dsn] = control(ds,ds0,Counter,Step_size_information,Do,x_plus,Path,Opt)
    if ~Do.stepback
        if ~Do.deflate
            if Counter.error == 0
                if Opt.step_size_control.iterations
                    solver_output = Step_size_information.current.solver_output;
                    dsn = step_size.iterations(ds,ds0,Counter,solver_output,Do,Path,Opt);
                elseif Opt.step_size_control.angle
                    solver_output = Step_size_information.current.solver_output;
                    dsn = step_size.angle(ds,ds0,Counter,solver_output,Do,Path,Opt);
                elseif Opt.step_size_control.curvature
                    if length(Path.l_all) > 3
                        solver_output = Step_size_information.current.solver_output;
                        dsn = step_size.curvature(ds,ds0,Counter,solver_output,Do,Path,Opt);
                    else
                        solver_output = Step_size_information.current.solver_output;
                        dsn = step_size.iterations(ds,ds0,Counter,solver_output,Do,Path,Opt);
                    end
                elseif Opt.step_size_control.pid
                    if length(Path.l_all) > 4
                        solver_output = Step_size_information.current.solver_output;
                        dsn = step_size.pid(ds,ds0,Counter,solver_output,Do,Path,Opt);
                    else
                        solver_output = Step_size_information.current.solver_output;
                        dsn = step_size.iterations(ds,ds0,Counter,solver_output,Do,Path,Opt);
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
        %% Limit to max./min. step size:
        dsn = min([norm(Opt.ds_max),dsn]);
        dsn = max([Opt.ds_min,dsn]);
    else
        xe = [Path.var_all(:,end);Path.l_all(end)];
        dsn = norm(x_plus-xe)/2;
    end
end