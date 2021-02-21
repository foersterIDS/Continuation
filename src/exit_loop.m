%% path continuation - exit_loop
%
%   Institute of Dynamics and Vibration Research
%   Leibniz University Hannover
%   03.11.2020 - Tido Kubatschek
%   21.02.2021 - Alwin Förster
%
function [do_continuation, exitflag, var_all, l_all, s_all, Opt] = exit_loop(do_continuation, exitflag, l_start, l_end, var_all, l_all, s_all, Opt, loop_counter, error_counter, bif_flag, bif, ds, fun_solution, solver_jacobian)
    %% eval. break function:
    try
        bfun = Opt.break_function(fun_solution,solver_jacobian,var_all(:,end),l_all(end));
    catch
        warning('Unable to evaluate user defined break function.');
        bfun = false;
    end
    %% exit without complete results:
    if error_counter>=Opt.max_error_counter
        exitflag = -1;
        warning('No valid result could be found for the last %d attempts.',Opt.max_error_counter);
        do_continuation = false;
    end
    %
    %% exit with l<l_start:
    if sign(l_end-l_start)*(l_all(end)-l_start)<0
        do_continuation = false;
        exitflag = 0;
    end
    %
    %% exit with success:
    if sign(l_end-l_start)*(l_all(end)-l_end)>=0
        do_continuation = false;
        exitflag = 1;
    end
    %
    %% exit with n_step_max reached:
    if loop_counter>=Opt.n_step_max
        do_continuation = false;
        exitflag = 2;
    end
    %% exit with bifurcation:
    if bif_flag>0 && Opt.stop_on_bifurcation
        do_continuation = false;
        exitflag = 3;
        var_all = var_all(:,1:bif(1,end));
        l_all = l_all(1:bif(1,end));
        s_all = s_all(1:bif(1,end));
    end
    %
    %% exit on closed curve:
    if Opt.closed_curve_detection
        [is_closed, Opt] = closed_curve(Opt,var_all,l_all,s_all, ds);
        if is_closed
            do_continuation = false;
            exitflag = 4;
            warning('closed curve detected. stopping continuation!');
        end
    end
    %
    %% exit due to break function:
    if bfun
        do_continuation = false;
        exitflag = 5;
    end
    %
end

