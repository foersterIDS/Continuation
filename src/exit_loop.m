%% path continuation - exit_loop
%
%   Institute of Dynamics and Vibration Research
%   Leibniz University Hannover
%   03.11.2020 - Tido Kubatschek
%   21.02.2021 - Alwin Förster
%
function [Do, Info, Path, break_fun_out, Opt] = exit_loop(Do, Info, l_start, l_end, Path, Opt, Counter, Bifurcation, ds, fun_solution, solver_jacobian, break_fun_out, val)
    %% eval. break function:
    %
    try
        if val
            [bfun,break_fun_out] = Opt.break_function(fun_solution,solver_jacobian,Path.var_all(:,end),Path.l_all(end),break_fun_out);
        else
            bfun = false;
        end
    catch
        warning('Unable to evaluate user defined break function.');
        bfun = false;
    end
    %
    %% exit with max_remove_counter reached:
    %
    if Counter.remove>Opt.max_remove_counter
        Do.continuation = false;
        Info.exitflag = -2;
    end
    %
    %% exit without complete results:
    %
    if Counter.error>=Opt.max_error_counter
        Info.exitflag = -1;
        warning('No valid result could be found for the last %d attempts.',Opt.max_error_counter);
        Do.continuation = false;
    end
    %
    %% exit with l<l_start:
    %
    if sign(l_end-l_start)*(Path.l_all(end)-l_start)<0
        Do.continuation = false;
        Info.exitflag = 0;
    end
    %
    %% exit with success:
    %
    if sign(l_end-l_start)*(Path.l_all(end)-l_end)>=0
        Do.continuation = false;
        Info.exitflag = 1;
    end
    %
    %% exit with n_step_max reached:
    %
    if Counter.loop>=Opt.n_step_max
        Do.continuation = false;
        Info.exitflag = 2;
    end
    %
    %% exit with bifurcation:
    %
    if Bifurcation.flag>0 && Opt.stop_on_bifurcation
        Do.continuation = false;
        Info.exitflag = 3;
        Path.var_all = Path.var_all(:,1:Bifurcation.bif(1,end));
        Path.l_all = Path.l_all(1:Bifurcation.bif(1,end));
        Path.s_all = Path.s_all(1:Bifurcation.bif(1,end));
    end
    %
    %% exit on closed curve:
    %
    if Opt.closed_curve_detection
        [is_closed, Opt] = closed_curve(Opt,Path, ds);
        if is_closed
            Do.continuation = false;
            Info.exitflag = 4;
            warning('closed curve detected. stopping continuation!');
        end
    end
    %
    %% exit due to break function:
    %
    if bfun
        Do.continuation = false;
        Info.exitflag = 5;
    end
    %
    %% exit due to max_continuation_steps reached
    %
    if Opt.max_continuation_steps<=numel(Path.l_all)
        Do.continuation = false;
        Info.exitflag = 6;
    end
    %
end

