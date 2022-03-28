%% path continuation - aux.exit_loop
%
%   Institute of Dynamics and Vibration Research
%   Leibniz University Hannover
%   03.11.2020 - Tido Kubatschek
%   21.02.2021 - Alwin Förster
%
function [Do, Info, Path, break_fun_out, Opt,Counter] = exit_loop(Do, Info, Path, Opt, Counter, Bifurcation, ds, fun_solution, Jacobian, break_fun_out, val)
    %% eval. break function:
    %
    try
        if val
            [bfun,break_fun_out] = Opt.break_function(fun_solution,Jacobian.solver,Path.var_all(:,end),Path.l_all(end),break_fun_out);
        else
            bfun = false;
        end
    catch
        aux.print_line(Opt,'--> Unable to evaluate user defined break function.\n');
        bfun = false;
    end
    %
    %% exit with max_remove_counter reached:
    %
    if Counter.remove>Opt.max_remove_counter
        Do.continuation = false;
        Info.exitflag = -3;
        Info.exit_msg = '--> continuation stoped: max_remove_counter has been reached.';
    end
    %
    %% exit without complete results:
    %
    if Counter.catch >= 3
        Info.exitflag = -2;
        Do.continuation = false;
        Info.exit_msg = sprintf('--> continuation stoped: Too many erros in function evaluation.');
    end
    %
    %% exit without complete results:
    %
    if Counter.error>=Opt.max_error_counter
        Info.exitflag = -1;
        Do.continuation = false;
        Info.exit_msg = sprintf('--> continuation stoped: No valid result could be found for the last %d attempts.',Opt.max_error_counter);
    end
    %
    %% exit with l<l_start:
    %
    if sign(Info.l_end-Info.l_start)*(Path.l_all(end)-Info.l_start)<0
        Do.continuation = false;
        Info.exitflag = 0;
        Info.exit_msg = '--> continuation stoped: l<l_start';
    end
    %
    %% exit with success:
    %
    if sign(Info.l_end-Info.l_start)*(Path.l_all(end)-Info.l_end)>=0 || sign(Opt.l_target-Info.l_start)*(Path.l_all(end)-Opt.l_target)>=0
        Do.continuation = false;
        Info.exitflag = 1;
        if sign(Info.l_end-Info.l_start)*(Path.l_all(end)-Info.l_end)>=0
            Info.exit_msg = '--> continuation completed: l_end reached';
        else
            Info.exit_msg = '--> continuation completed: l_target reached';
        end
    end
    %
    %% exit with n_step_max reached:
    %
    if Counter.loop>=Opt.n_step_max
        Do.continuation = false;
        Info.exitflag = 2;
        Info.exit_msg = '--> continuation completed: n_step_max reached';
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
        Info.exit_msg = '--> continuation completed: bifurcation reached';
    end
    %
    %% exit on closed curve:
    %
    if Opt.closed_curve_detection
        [is_closed, Opt, Counter] = aux.closed_curve(Opt,Path,ds,Counter);
        if is_closed
            Do.continuation = false;
            Info.exitflag = 4;
            Info.exit_msg = '--> continuation completed: closed curve detected';
        end
    end
    %
    %% exit due to break function:
    %
    if bfun
        Do.continuation = false;
        Info.exitflag = 5;
        Info.exit_msg = '--> continuation completed: user-defined break function';
    end
    %
    %% exit due to user input:
    %
    if Do.stop_manually
        Do.continuation = false;
        Info.exitflag = 6;
        Info.exit_msg = '--> continuation stoped by user';
    end
    %
end

