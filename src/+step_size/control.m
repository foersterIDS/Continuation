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
%   08.05.2020 - Alwin Förster
%   03.06.2020 - Niklas Marhenke
%   21.10.2020 - Tido Kubatschek
%
function [dsn,Counter,event_out,Opt] = control(ds,Counter,Solver,Do,Plus,Path,Jacobian,Opt,Info,event_out,Initial)
    if ~Do.stepback
        if ~Do.deflate
            if Counter.error == 0
                if Opt.step_size_event
                    [dsn,Counter,event_out,changed,Opt] = step_size.event_adjustment(ds,Path,Counter,Opt,event_out,Initial);
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
                            xi = step_size.angle_change(Solver,Path,Opt);
                        else
                            xi = step_size.iterations_polynomial(Solver,Opt);
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
                            xi = step_size.angle_custom(Solver,Path,Opt);
                        else
                            xi = step_size.iterations_polynomial(Solver,Opt);
                        end
                    %
                    % contraction
                    %
                    elseif Opt.step_size_control.contraction
                        xi = step_size.contraction(Solver,Opt);
                    %
                    % error
                    %
                    elseif Opt.step_size_control.error
                        xi = step_size.error(Solver,Path,Opt);
                    %
                    % error_alt
                    %
                    elseif Opt.step_size_control.error_alt
                        xi = step_size.error_alt(Solver,Path,Opt);
                    %
                    % fayezioghani
                    %
                    elseif Opt.step_size_control.fayezioghani
                        %
                        % Check if there are enough solution points to use
                        % method. Otherwise use iterations.
                        %
                        if length(Path.l_all) > 2
                            xi = step_size.fayezioghani(ds,Solver,Path,Jacobian.last,Opt);
                        else
                            xi = step_size.iterations_polynomial(Solver,Opt);
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
                        xi = step_size.iterations_exponential(Solver,Opt);
                    % 
                    % iterations based - polynomial
                    %
                    elseif Opt.step_size_control.iterations_polynomial
                        xi = step_size.iterations_polynomial(Solver,Opt);
                    % 
                    % multiplicative method
                    %
                    elseif Opt.step_size_control.multiplicative
                        xi = step_size.multiplicative(Solver,Path,Opt);
                    % 
                    % multiplicative_alt method
                    %
                    elseif Opt.step_size_control.multiplicative_alt
                        xi = step_size.multiplicative_alt(Solver,Path,Opt);
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
                            xi = step_size.iterations_polynomial(Solver,Opt);
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
                            xi = step_size.iterations_polynomial(Solver,Opt);
                        end
                    %
                    % szyszkowski
                    %
                    elseif Opt.step_size_control.szyszkowski
                        if length(Path.l_all) > 3
                            xi = step_size.szyszkowski(Solver,Path,Opt);
                        else
                            xi = step_size.iterations_polynomial(Solver,Opt);
                        end
                    %
                    % yoon
                    %
                    elseif Opt.step_size_control.yoon
                        if length(Path.l_all) > 3
                            xi = step_size.yoon(Solver,Path,Opt);
                        else
                            xi = step_size.iterations_polynomial(Solver,Opt);
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
                    %% Limit to max./min. step size, also limit to Opt.max_step_size_change*ds / ds/Opt.max_step_size_change*ds:
                    dsn = min([norm(Opt.ds_max),Opt.max_step_size_change*ds,dsn]);
                    dsn = max([Opt.ds_min,ds/Opt.max_step_size_change,dsn]);
                end
            else
                if Opt.step_size_control.fix
                    dsn = Info.ds0;
                else
                    dsn = ds/2;
                    %% Limit to max./min. step size, also limit to Opt.max_step_size_change*ds / ds/Opt.max_step_size_change*ds:
                    dsn = min([norm(Opt.ds_max),Opt.max_step_size_change*ds,dsn]);
                    dsn = max([Opt.ds_min,ds/Opt.max_step_size_change,dsn]);
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
            dsn = norm(Plus.x-xe)/2;
        end
    end
end