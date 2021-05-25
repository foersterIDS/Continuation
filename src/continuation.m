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
function [var_all,l_all,exitflag,bif,s_all,last_jacobian,break_fun_out] = continuation(fun,var0,l_start,l_end,ds0,varargin)
    %% initialize
    %
    exitflag = -1;
    warning on;
    Opt = continuation_input(varargin,fun,var0,l_start,l_end);
    [solver,predictor_solver,default_solver_output] = continuation_solver(Opt);
    res_arle = residual_corrector(Opt);
    ds = ds0;
    nv = length(var0);
    catch_counter = 0;
    catch_counter_old = 0;
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
    [var_all,fun_initial,initial_exitflag,solver_output,initial_jacobian] = solver(residual_initial,var0,Opt.dscale0(1:end-1));
    solver_jacobian = initial_jacobian;
    bif = [];
    x_plus = [];
    last_jacobian = [];
    break_fun_out = [];
    if initial_exitflag>0
        l_all = Opt.l_0;
        s_all = 0;
        do_continuation = true;
        do_loop = true;
        solver_jacobian = [solver_jacobian,numeric_jacobian(@(x) fun(x(1:nv),x(nv+1)),[var_all;Opt.l_0],'central_value',fun_initial,'derivative_dimensions',nv+1,'diffquot',Opt.diffquot)];
        previous_jacobian = solver_jacobian;
        last_jacobian = solver_jacobian;
        [~,break_fun_out] = Opt.break_function(fun_initial,solver_jacobian,var_all,l_all,break_fun_out);
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
        catch_counter_old = catch_counter;
        %
        %% residual
        %
        residual = @(x) merge_residuals(fun,res_arle,x,[var_all;l_all],ds,Opt);
        if do_deflate
            try
                residual = @(x) deflation(residual,x_deflation,x,Opt);
            catch
                error('Error occured during deflation.');
                catch_counter = catch_counter + 1;
                if catch_counter >= 3
                    warning('Error in input! catch was used too often!');
                    break;
                end
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
                [var_predictor,l_predictor,fun_predictor,s_predictor] = predictor(var_all,l_all,s_all,ds,last_jacobian,fun,res_arle,predictor_solver,Opt);
                x_predictor = [var_predictor;l_predictor];
            end
        catch
            [var_predictor,l_predictor,fun_predictor,s_predictor] = predictor(var_all,l_all,s_all,ds,last_jacobian,fun,res_arle,predictor_solver,Opt);
            x_predictor = [var_predictor;l_predictor];
            warning('predictor: catch!');
            catch_counter = catch_counter + 1;
            if catch_counter >= 3
                    warning('Error in input! catch was used too often!');
                    break;
            end
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
                %            
                if Opt.solver.fsolve && Opt.solver_force1it
                    [x_solution,fun_solution,solver_exitflag,solver_output,solver_jacobian] = solver(residual,x_predictor,dscale);
                    if solver_output.iterations < 1
                        % perturbate initial solution by tolerance of
                        % solver
                        pert = Opt.solver_tol * ones(numel(x_predictor),1) / numel(x_predictor);
                        x_predictor = x_predictor + pert;
                        [x_solution,fun_solution,solver_exitflag,solver_output,solver_jacobian] = solver(residual,x_predictor,dscale);
                    end
                else
                    [x_solution,fun_solution,solver_exitflag,solver_output,solver_jacobian] = solver(residual,x_predictor,dscale);
                end
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
            catch_counter = catch_counter + 1;
            if catch_counter >= 3
                    warning('Error in input! catch was used too often!');
                    break;
            end
        end
        %
        %% check result
        %
        [val,is_reverse, catch_flag] = validate_result(x_solution,fun_solution,var_all,l_all,ds,solver_exitflag,solver_jacobian,last_jacobian,do_convergeToTarget,Opt);
        if val
            %% valid result
            if isempty(x_plus)
                var_all = [var_all,x_solution(1:end-1)];
                l_all = [l_all,x_solution(end)];
                s_all = [s_all,s_all(end)+norm(x_solution-[var_all(:,end-1);l_all(end-1)])];
                previous_jacobian = last_jacobian;
                last_jacobian = solver_jacobian;
            else
                var_all = [var_all,x_solution(1:end-1),x_plus(1:end-1)];
                l_all = [l_all,x_solution(end),x_plus(end)];
                s_all = [s_all,s_all(end)+norm(x_solution-[var_all(:,end-2);l_all(end-2)])*[1,1]+norm(x_plus-x_solution)*[0,1]];
                x_plus = [];
                previous_jacobian = solver_jacobian;
                last_jacobian = plus_jacobian;
                plus_jacobian = [];
            end
            do_deflate = false;
            do_homotopy = false;
            do_stepback = false;
            error_counter = 0;
            step_loop = step_loop + 1;
            if do_convergeToTarget
                last_jacobian = get_jacobian(fun,var_all(:,end),l_all(end),Opt);
            end
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
                plus_jacobian = last_jacobian;
                last_jacobian = previous_jacobian;
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
            if catch_flag
                catch_counter = catch_counter + 1;
                if catch_counter >= 3
                    warning('Error in input! catch was used too often!');
                    break;
                end
            end
        end
        %
        %% Bifurcations
        %
        bif_flag = 0;
        if ison(Opt.bifurcation) && val && ~do_homotopy && numel(l_all)>2
            if ~is_current_jacobian
                %% get jacobian if not current
                solver_jacobian = get_jacobian(fun,var_all(:,end),l_all(end),Opt);
            end
            [bif,sign_det_jacobian,bif_flag,var_all,l_all,s_all] = check_bifurcation(fun,solver_jacobian(1:nv,1:nv),var_all,l_all,s_all,bif,sign_det_jacobian,res_arle,predictor_solver,Opt);
        elseif ison(Opt.bifurcation) && val && numel(l_all)<=2
            if ~is_current_jacobian
                %% get jacobian if not current
                solver_jacobian = get_jacobian(fun,var_all(:,end),l_all(end),Opt);
            end
            sign_det_jacobian = sign(det(solver_jacobian(1:nv,1:nv)));
        end
        %
        %% step size control
        %
        dsim1 = ds;
        ds = step_size_control(ds,ds0,error_counter,solver_output,do_deflate,do_stepback,x_plus,var_all,l_all,s_all,Opt);
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
        [do_continuation,exitflag,var_all,l_all,s_all,break_fun_out,Opt] = exit_loop(do_continuation,exitflag,l_start,l_end,var_all,l_all,s_all,Opt,loop_counter,error_counter,bif_flag,bif,ds,fun_solution,solver_jacobian,break_fun_out,val);
        %
        %% live plot
        %
        if ison(Opt.plot) && val
            try
                [pl, Opt] = live_plot(Opt, nv, l_start, l_end, l_all, var_all, s_all, ds, dsim1, solver_output.iterations, loop_counter, fun_predictor, s_predictor, pl, bif_flag, bif);
            catch
                warning('The plot update has failed.');
            end
        end
        %
    end
    %
    %% live plot finalization
    %
    if ison(Opt.plot) && initial_exitflag>0
        try
            live_plot(Opt, nv, l_start, l_end, l_all, var_all, s_all, ds, dsim1, solver_output.iterations, loop_counter, fun_predictor, s_predictor, pl, -1);
        catch
            warning('The plot update has failed.');
        end
    end
    %
    %% bifurcation tracing
    %
    if Opt.bifurcation.trace
        try
            [var_all,l_all,s_all,bif] = trace_bifurcations(Opt,var_all,l_all,s_all,bif,solver,fun,l_start,l_end,res_arle,predictor_solver);
            last_jacobian = [];
        catch
            warning('Failed to trace bifurcations.');
        end
    end
    %
    %% final disp
    %
    if Opt.display
        fprintf('Time Elapsed: %.3f s\n',toc(t_display));
    end
    %
    %% catch counter
    if catch_counter_old == catch_counter
        catch_counter = 0;
    end
end