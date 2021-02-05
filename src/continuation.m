%% path continuation - continuation
%
%   Institute of Dynamics and Vibration Research
%   Leibniz University Hannover
%   08.05.2020 - Alwin Förster
%
%   [var_all,l_all,exitflag,bif] = continuation(fun,var0,l_start,l_end,ds0,varargin)
%
%   fun = fun(var,l) != 0
%   l_start <= l <= l_end
%   ds0: initial stepsize
%
function [var_all,l_all,exitflag,bif,s_all] = continuation(fun,var0,l_start,l_end,ds0,varargin)
    %% initialize
    %
    exitflag = -1;
    warning on;
    Opt = continuation_input(varargin,fun,var0,l_start,l_end);
    [solver,predictor_solver,default_solver_output] = continuation_solver(Opt);
    res_arle = residual_arclength(Opt);
    ds = ds0;
    nv = length(var0);
    do_deflate = false;
    do_homotopy = false;
    do_stepback = false;
    do_convergeToTarget = false;
    error_counter = 0;
    loop_counter = 0;
    step_loop = 0;
    if Opt.display
        fprintf('Starting path continuation...\n');
        t_display = tic;
    end
    %
    %% find initial solution
    %
    residual_initial = @(v) residual_fixed_value(fun,v,Opt.l_0,Opt);
    [var_all,~,initial_exitflag,solver_output,initial_jacobian] = solver(residual_initial,var0,Opt.dscale0(1:end-1));
    solver_jacobian = initial_jacobian;
    bif = [];
    x_plus = [];
    if initial_exitflag>0
        l_all = Opt.l_0;
        s_all = 0;
        do_continuation = true;
        do_loop = true;
        if Opt.display
            fprintf('Initial solution at l = %.2e\n',Opt.l_0);
        end
        if ison(Opt.bifurcation)
            sign_det_jacobian = sign(det(initial_jacobian));
        end
        if ison(Opt.plot)
            [pl, Opt] = live_plot(Opt, nv, l_start, l_end, l_all, var_all, s_all, ds0, ds0, solver_output.iterations, loop_counter);
        end
    else
        var_all = [];
        l_all = [];
        s_all = [];
        exitflag = -2; % no initial solution found
        do_continuation = false;
        if Opt.display
            fprintf('No initial solution found.\n');
        end
    end
    %
    %% continuation
    %
    while do_continuation
        %% initialize loop
        %
        loop_counter = loop_counter+1;
        is_current_jacobian = false;
        %
        %% residual
        %
        residual = @(x) merge_residuals(fun,res_arle,x,[var_all;l_all],ds,Opt);
        if do_deflate
            try
                residual = @(x) deflation(residual,x_deflation,x,Opt);
            catch
                error('Error occured during deflation.');
            end
        end
        %
        %% predictor
        %
        try
            if do_homotopy
                %% Homotopy
                x_last_step = [var_all(:,end);l_all(end)];
                x_predictor = homotopy(residual,x_last_step,Opt);
                var_predictor = x_predictor(1:end-1);
                l_predictor = x_predictor(end);
            else
                [var_predictor,l_predictor,fun_predictor,s_predictor] = predictor(var_all,l_all,s_all,ds,solver_jacobian,fun,res_arle,predictor_solver,Opt);
                x_predictor = [var_predictor;l_predictor];
            end
        catch
            [var_predictor,l_predictor,fun_predictor,s_predictor] = predictor(var_all,l_all,s_all,ds,solver_jacobian,fun,res_arle,predictor_solver,Opt);
            x_predictor = [var_predictor;l_predictor];
            warning('predictor: catch!');
        end
        %
        %% solve
        %
        try
            dscale = get_dscale(Opt,var_all,l_all);
            if sign(l_all(end)-Opt.l_target)*sign(l_predictor-Opt.l_target)<=0
                %% try to converge to target
                residual_target = @(v) residual_fixed_value(fun,v,Opt.l_target,Opt);
                var_predictor_ctt = (var_predictor - var_all(:,end))*(abs(Opt.l_target-l_all(end))/abs(l_predictor-l_all(end)))+var_all(:,end);
                [var_solution,fun_solution,solver_exitflag,solver_output,solver_jacobian] = solver(residual_target,var_predictor_ctt,dscale(1:end-1));
                x_solution = [var_solution;Opt.l_target];
                do_convergeToTarget = true;
            else
                %% regular solver
                [x_solution,fun_solution,solver_exitflag,solver_output,solver_jacobian] = solver(residual,x_predictor,dscale);
                do_convergeToTarget = false;
            end
            is_current_jacobian = true;
        catch
            x_solution = NaN(size(x_predictor));
            fun_solution = inf(size(x_predictor));
            solver_exitflag = -2;
            solver_output = default_solver_output;
            do_convergeToTarget = false;
            warning('solve: catch!');
        end
        %
        %% check result
        %
        [val,is_reverse] = validate_result(x_solution,fun_solution,var_all,l_all,ds,solver_exitflag,do_convergeToTarget,Opt);
        if val
            %% valid result
            if isempty(x_plus)
                var_all = [var_all,x_solution(1:end-1)];
                l_all = [l_all,x_solution(end)];
                s_all = [s_all,s_all(end)+norm(x_solution-[var_all(:,end-1);l_all(end-1)])];
            else
                var_all = [var_all,x_solution(1:end-1),x_plus(1:end-1)];
                l_all = [l_all,x_solution(end),x_plus(end)];
                s_all = [s_all,s_all(end)+norm(x_solution-[var_all(:,end-2);l_all(end-2)])*[1,1]+norm(x_plus-x_solution)*[0,1]];
                x_plus = [];
            end
            do_deflate = false;
            do_homotopy = false;
            do_stepback = false;
            error_counter = 0;
            step_loop = step_loop + 1;
        else
            %% invalid result
            error_counter = error_counter+1;
            if Opt.deflation && ~do_deflate && ~isnan(sum(x_solution(:,end))) && is_reverse && error_counter>=Opt.deflation_error_counter
                do_deflate = true;
                x_deflation = x_solution;
            else
                do_deflate = false;
            end
            do_stepback = (error_counter==Opt.stepback_error_counter) && (length(l_all)>1);
            if (error_counter==Opt.stepback_error_counter) && (length(l_all)>1)
                x_plus = [var_all(:,end);l_all(end)];
                var_all(:,end) = [];
                l_all(end) = [];
                s_all(end) = [];
            elseif (error_counter==Opt.stepback_error_counter+1) && (length(l_all)>1)
                var_all = [var_all,x_plus(1:end-1)];
                l_all = [l_all,x_plus(end)];
                s_all = [s_all,s_all(end)+norm([var_all(:,end);l_all(end)]-[var_all(:,end-1);l_all(end-1)])];
                x_plus = [];
            else
                x_plus = [];
            end
            if Opt.include_reverse && is_reverse && solver_exitflag>0
                [var_all,l_all] = include_reverse(x_solution,var_all,l_all);
            end
            if ison(Opt.homotopy) && error_counter>=Opt.homotopy_error_counter
                do_homotopy = true;
            else
                do_homotopy = false;
            end
        end
        %
        %% Bifurcations
        %
        bif_flag = 0;
        if ison(Opt.bifurcation) && val && ~do_homotopy && numel(l_all)>2
            if ~is_current_jacobian
                %% get jacobian if not current
                solver_jacobian = get_jacobian(fun,var_all(:,end),l_all(end));
            end
            [bif,sign_det_jacobian,bif_flag,var_all,l_all,s_all] = check_bifurcation(fun,solver_jacobian(1:nv,1:nv),var_all,l_all,s_all,bif,sign_det_jacobian,res_arle,predictor_solver,Opt);
        elseif ison(Opt.bifurcation) && val && numel(l_all)<=2
            if ~is_current_jacobian
                %% get jacobian if not current
                solver_jacobian = get_jacobian(fun,var_all(:,end),l_all(end));
            end
            sign_det_jacobian = sign(det(solver_jacobian(1:nv,1:nv)));
        end
        %
        %% step size control
        %
        dsim1 = ds;
        ds = step_size_control(ds,ds0,error_counter,solver_output,do_deflate,do_stepback,x_plus,var_all,l_all,s_all,Opt);
        %
        %% live plot
        %
        if ison(Opt.plot) && val
            [pl, Opt] = live_plot(Opt, nv, l_start, l_end, l_all, var_all, s_all, ds, dsim1, solver_output.iterations, loop_counter, fun_predictor, s_predictor, pl, bif_flag, bif);
        end
        %
        %% end loop
        %
        if Opt.display
            if val
                fprintf('-----> continued at l = %.4e\t|\tnew arc-length: ds = %.2e\t|\tloop counter = %d\t|\tstep = %d\t|\titerations = %d/%d\n',l_all(end),ds,loop_counter,step_loop,solver_output.iterations,Opt.n_iter_opt);
            else
                fprintf('-----> invalid point\t\t\t\t|\tnew arc-length: ds = %.2e\t|\tloop counter = %d\t|\tstep = %d\t|\titerations = %d/%d\n',ds,loop_counter,step_loop,solver_output.iterations,Opt.n_iter_opt);
            end
        end
        [do_continuation, exitflag, var_all, l_all, s_all, Opt] = exit_loop(do_continuation, exitflag, l_start, l_end, var_all, l_all, s_all, Opt, loop_counter, error_counter, bif_flag, bif, ds);
        %
    end
    %
    %% live plot finalization
    %
    if ison(Opt.plot) && initial_exitflag>0
        live_plot(Opt, nv, l_start, l_end, l_all, var_all, s_all, ds, dsim1, solver_output.iterations, loop_counter, fun_predictor, s_predictor, pl, -1);
    end
    %
    %% bifurcation tracing
    %
    if Opt.bifurcation.trace
        [var_all,l_all,s_all,bif] = trace_bifurcations(Opt,var_all,l_all,s_all,bif,solver,fun,l_start,l_end,res_arle,predictor_solver);
    end
    %
    %% final disp
    %
    if Opt.display
        fprintf('Time Elapsed: %.3f s\n',toc(t_display));
    end
    %
end