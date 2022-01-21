%% path continuation - step_size.control
%
%   Institute of Dynamics and Vibration Research
%   Leibniz University Hannover
%   08.05.2020 - Alwin Förster
%   03.06.2020 - Niklas Marhenke
%   21.10.2020 - Tido Kubatschek
%
function [dsn] = control(ds,Counter,solver_output,Do,x_plus,Path,last_jacobian,Opt,Info)
    if ~Do.stepback
        if ~Do.deflate
            if Counter.error == 0
                %% Determine step size control method
                %
                %
                % angle change based
                %
                if Opt.step_size_control.angle_change
                    %
                    % Check if there are enough solution points to use
                    % method. Otherwise use iterations.
                    %
                    if length(Path.l_all) > 3
                        xi = step_size.angle_change(solver_output,Path,Opt);
                    else
                        xi = step_size.iterations_polynomial(solver_output,Opt);
                    end
                %
                % angle custom
                %
                elseif Opt.step_size_control.angle_custom
                    %
                    % Check if there are enough solution points to use
                    % method. Otherwise use iterations.
                    %
                    if length(Path.l_all) > 2
                        xi = step_size.angle_custom(solver_output,Path,Opt);
                    else
                        xi = step_size.iterations_polynomial(solver_output,Opt);
                    end
                %
                % contraction
                %
                elseif Opt.step_size_control.contraction
                    xi = step_size.contraction(solver_output);
                %
                % error
                %
                elseif Opt.step_size_control.error
                    xi = step_size.error(solver_output,Path,Opt);
                %
                % fayezioghani
                %
                elseif Opt.step_size_control.fayezioghani
                    %
                    % Check if there are enough solution points to use
                    % method. Otherwise use iterations.
                    %
                    if length(Path.l_all) > 2
                        xi = step_size.fayezioghani(ds,solver_output,Path,last_jacobian,Opt);
                    else
                        xi = step_size.iterations_polynomial(solver_output,Opt);
                    end
                %
                % fixed step size
                %
                elseif Opt.step_size_control.fix
                    xi = 1;
                % 
                % iterations based - exponential
                %
                elseif Opt.step_size_control.iterations_exponential
                    xi = step_size.iterations_exponential(solver_output,Opt);
                % 
                % iterations based - polynomial
                %
                elseif Opt.step_size_control.iterations_polynomial
                    xi = step_size.iterations_polynomial(solver_output,Opt);
                % 
                % multiplicative method
                %
                elseif Opt.step_size_control.multiplicative
                    xi = step_size.multiplicative(solver_output,Path,Opt);
                %
                % pid control - custom
                %
                elseif Opt.step_size_control.pid_custom
                    %
                    % Check if there are enough solution points to use
                    % method. Otherwise use iterations.
                    %
                    if length(Path.l_all) > 4
                        xi = step_size.pid_custom(Path,Opt);
                    else
                        xi = step_size.iterations_polynomial(solver_output,Opt);
                    end
                %
                % pid control - valli
                %
                elseif Opt.step_size_control.pid_valli
                    %
                    % Check if there are enough solution points to use
                    % method. Otherwise use iterations.
                    %
                    if length(Path.l_all) > 4
                        xi = step_size.pid_valli(Path,Opt);
                    else
                        xi = step_size.iterations_polynomial(solver_output,Opt);
                    end
                %
                % szyszkowski
                %
                elseif Opt.step_size_control.szyszkowski
                    if length(Path.l_all) > 3
                        xi = step_size.szyszkowski(solver_output,Path,Opt);
                    else
                        xi = step_size.iterations_polynomial(solver_output,Opt);
                    end
                %
                % yoon
                %
                elseif Opt.step_size_control.yoon
                    if length(Path.l_all) > 3
                        xi = step_size.yoon(solver_output,Path,Opt);
                    else
                        xi = step_size.iterations_polynomial(solver_output,Opt);
                    end
                %
                % wrong method
                %
                else
                    error('Invalid settings for step size control!');
                end
                %
                %% adapt step size
                dsn = xi * ds;
            else
                if Opt.step_size_control.fix
                    dsn = Info.ds0;
                else
                    dsn = ds/2;
                end
            end
        else
            dsn = ds;
        end
        %% Limit to max./min. step size, also limit to 2*ds/0.5*ds:
        dsn = min([norm(Opt.ds_max),2*ds,dsn]);
        dsn = max([Opt.ds_min,ds/2,dsn]);
    else
        if Opt.step_size_control.fix
            dsn = Info.ds0;
        else
            xe = [Path.var_all(:,end);Path.l_all(end)];
            dsn = norm(x_plus-xe)/2;
        end
    end
end