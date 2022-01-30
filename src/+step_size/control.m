%% path continuation - step_size.control
%  Adjusts stepsize by choosing an adaption method. 
%
%
%
%  See the <a href="matlab:open('..\doc\html\continuation.html')">documentation</a> or
%  see the control methods:
%  -- <a href="matlab:doc('step_size.angle_change')">step_size.angle_change</a>
%  -- <a href="matlab:doc('step_size.angle_custom')">step_size.angle_custom</a>
%  -- <a href="matlab:doc('step_size.contraction')">step_size.contraction</a>
%  -- <a href="matlab:doc('step_size.error')">step_size.error</a>
%  -- <a href="matlab:doc('step_size.event_adjustment')">step_size.event_adjustment</a>
%  -- <a href="matlab:doc('step_size.fayezioghani')">step_size.fayezioghani</a>
%  -- <a href="matlab:doc('step_size.iterations_exponential')">step_size.iterations_exponential</a>
%  -- <a href="matlab:doc('step_size.iterations_polynomial')">step_size.iterations_polynomial</a>
%  -- <a href="matlab:doc('step_size.multiplicative')">step_size.multiplicative</a>
%  -- <a href="matlab:doc('step_size.pid_custom')">step_size.pid_custom</a>
%  -- <a href="matlab:doc('step_size.pid_valli')">step_size.pid_valli</a>
%  -- <a href="matlab:doc('step_size.szyszkowski')">step_size.szyszkowski</a>
%  -- <a href="matlab:doc('step_size.yoon')">step_size.yoon</a>
%
%
%   Institute of Dynamics and Vibration Research
%   Leibniz University Hannover
%   08.05.2020 - Alwin F�rster
%   03.06.2020 - Niklas Marhenke
%   21.10.2020 - Tido Kubatschek
%
function [dsn,Counter,event_out] = control(ds,Counter,solver_output,Do,x_plus,Path,last_jacobian,Opt,Info,event_out)
    if ~Do.stepback
        if ~Do.deflate
            if Counter.error == 0
                if Opt.step_size_event
                    [dsn,Counter,event_out,changed] = step_size.event_adjustment(ds,Path,Counter,Opt,event_out);
                else
                    changed = false;
                end
                if ~changed
                    %
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
                        xi = step_size.contraction(solver_output,Opt);
                    %
                    % error
                    %
                    elseif Opt.step_size_control.error
                        xi = step_size.error(solver_output,Path,Opt);
                    %
                    % error_alt
                    %
                    elseif Opt.step_size_control.error_alt
                        xi = step_size.error_alt(solver_output,Path,Opt);
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
                    % multiplicative_alt method
                    %
                    elseif Opt.step_size_control.multiplicative_alt
                        xi = step_size.multiplicative_alt(solver_output,Path,Opt);
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
                    %% Limit to max./min. step size, also limit to 2*ds/0.5*ds:
                    dsn = min([norm(Opt.ds_max),2*ds,dsn]);
                    dsn = max([Opt.ds_min,ds/2,dsn]);
                end
            else
                if Opt.step_size_control.fix
                    dsn = Info.ds0;
                else
                    dsn = ds/2;
                    %% Limit to max./min. step size, also limit to 2*ds/0.5*ds:
                    dsn = min([norm(Opt.ds_max),dsn]);
                    dsn = max([Opt.ds_min,dsn]);
                end
            end
        else
            dsn = ds;
        end
    else
        if Opt.step_size_control.fix
            dsn = Info.ds0;
        else
            xe = [Path.var_all(:,end);Path.l_all(end)];
            dsn = norm(x_plus-xe)/2;
        end
    end
end