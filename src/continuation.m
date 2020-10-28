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
    Opt = continuation_input(varargin,fun,var0,l_start,l_end);
    [solver,default_solver_output] = continuation_solver(Opt);
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
    [var_all,~,initial_exitflag,~,initial_jacobian] = solver(residual_initial,var0);
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
        if Opt.plot
            pl = live_plot(Opt, nv, l_start, l_end, l_all, var_all);
        end
    else
        var_all = [];
        l_all = [];
        s_all = [];
        exitflag = -2; % no initial solution found
        do_continuation = false;
        do_loop = false;
        if Opt.display
            fprintf('No initial solution found.\n');
        end
    end
    %
    %% continuation
    %
    while do_continuation && ~Opt.unique
        %% initialize loop
        %
        loop_counter = loop_counter+1;
        is_current_jacobian = false;
        %
        %% residual and predictor
        %
        residual = @(x) merge_residuals(fun,res_arle,x,[var_all;l_all],ds,Opt);
        if do_deflate
            try
                residual = @(x) deflation(residual,x_deflation,x,Opt);
            catch
                error('Error occured during deflation.');
            end
        end
        try
            if do_homotopy
                %% Homotopy
                x_last_step = [var_all(:,end);l_all(end)];
                x_predictor = homotopy(residual,x_last_step,Opt);
            else
                [var_predictor,l_predictor] = predictor(var_all,l_all,s_all,ds,Opt);
                x_predictor = [var_predictor;l_predictor];
            end
        catch
            [var_predictor,l_predictor] = predictor(var_all,l_all,s_all,ds,Opt);
            x_predictor = [var_predictor;l_predictor];
        end
        %
        %% solve
        %
        try
            if sign(l_all(end)-Opt.l_target)*sign(l_predictor-Opt.l_target)<=0
                %% try to converge to target
                residual_target = @(v) residual_fixed_value(fun,v,Opt.l_target,Opt);
                [var_solution,fun_solution,solver_exitflag,solver_output,solver_jacobian] = solver(residual_target,var_predictor);
                x_solution = [var_solution;Opt.l_target];
                do_convergeToTarget = true;
            else
                %% regular solver
                [x_solution,fun_solution,solver_exitflag,solver_output,solver_jacobian] = solver(residual,x_predictor);
                do_convergeToTarget = false;
            end
            is_current_jacobian = true;
        catch
            x_solution = NaN(size(x_predictor));
            fun_solution = inf(size(x_predictor));
            solver_exitflag = -2;
            solver_output = default_solver_output;
            do_convergeToTarget = false;
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
            [bif,sign_det_jacobian,bif_flag,var_all,l_all,s_all] = check_bifurcation(fun,solver_jacobian(1:nv,1:nv),var_all,l_all,s_all,bif,sign_det_jacobian,Opt);
        elseif ison(Opt.bifurcation) && val && numel(l_all)<=2
            if ~is_current_jacobian
                %% get jacobian if not current
                solver_jacobian = get_jacobian(fun,var_all(:,end),l_all(end));
            end
            sign_det_jacobian = sign(det(solver_jacobian(1:nv,1:nv)));
        end
        %
        %% adjust arc-length
        %
        ds = step_size_control(ds,ds0,error_counter,solver_output,do_deflate,do_stepback,x_plus,var_all,l_all,s_all,Opt);
        %
        %% live plot
        %
        if Opt.plot && val
            pl = live_plot(Opt, nv, l_start, l_end, l_all, var_all, pl, bif_flag, bif);
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
        % exit with success:
        if sign(l_end-l_start)*(l_all(end)-l_end)>=0
            do_continuation = false;
            exitflag = 1;
        end
        % exit with n_step_max reached:
        if loop_counter>=Opt.n_step_max
            do_continuation = false;
            exitflag = 2;
        end
        % exit with l<l_start:
        if sign(l_end-l_start)*(l_all(end)-l_start)<0
            do_continuation = false;
            exitflag = 0;
        end
        % exit without complete results:
        if error_counter>=Opt.max_error_counter
            exitflag = -1;
            warning('No valid result could be found for the last %d attempts.',Opt.max_error_counter);
            do_continuation = false;
        end
        % exit with bifurcation:
        if bif_flag>0 && Opt.stop_on_bifurcation
            do_continuation = false;
            exitflag = 3;
            var_all = var_all(:,1:bif(1,end));
            l_all = l_all(1:bif(1,end));
            s_all = s_all(1:bif(1,end));
        end
        %
    end
    %
    %% bifurcation tracing
    %
    if Opt.bifurcation.trace && ~Opt.unique
        % TODO: start new continuation at bifurcations
        [var_all,l_all,s_all,bif] = trace_bifurcations(Opt,var_all,l_all,s_all,bif,solver,fun,l_start,l_end);
    end
    %
    %% unique - Durchlaufen mit festen Werten für l
    %
    if Opt.unique && do_loop
        do_loop = true;
        l_all = l_start;
        s_all = 0;
        l_v = l_start:ds0:l_end;
        %
        while loop_counter < length(l_v) && do_loop
            %% initialize loop
            %
            loop_counter = loop_counter + 1;
            l_aktuell = l_v(loop_counter);
            residual = @(v) residual_fixed_value(fun,v,l_aktuell,Opt);
            %
            %% find predictor
            %
            if sum(isnan(var_all(:,loop_counter))) > 0 % no last solution found
                var_aktuell = var0; % try with first solution
            else
                var_aktuell = var_all(:,loop_counter); % else try with last solution
            end
            %
            %% get next solution
            %
            [var_aktuell,~,exitflag,~,~] = solver(residual,var_aktuell);
            %
            %% check exitflag
            %
            if exitflag > 0
                step_loop = step_loop + 1;
                if Opt.display
                    fprintf('-----> Found solution at l = %.4e\t|\tloop counter = %d\t|\tstep = %d\n',l_aktuell, loop_counter, step_loop);
                end
            else
                var_aktuell = NaN(size(var_aktuell));
                error_counter = error_counter + 1;
                if Opt.display
                    fprintf('-----> No solution found at l = %.4e\t|\tloop counter = %d\t|\tstep = %d\n', l_aktuell, loop_counter, step_loop);
                end
            end
            %
            %% add solution 
            %
            var_all = [var_all, var_aktuell];
            l_all = [l_all, l_aktuell];
            s_all = [s_all,s_all(end)+norm([var_all(:,end);l_all(end)]-[var_all(:,end-1);l_all(end-1)])];
            %
            %% check for bifurcations
            %
% TODO:
            bif_flag = -2;
            bif = [];
            %
            %% check for error counter
            %
            if error_counter >= Opt.max_error_counter
                exitflag = -1;
                warning('No valid result could be found for the last %d attempts.',Opt.max_error_counter);
                do_loop = false;
            end
            %
            %% live plot
            %
            if Opt.plot
                pl = live_plot(Opt, nv, l_start, l_end, l_all, var_all, pl, bif_flag, bif);
            end
            %
        end
        %
    end
    %
    %% live plot finalization
    %
    if Opt.plot && initial_exitflag>0
        live_plot(Opt, nv, l_start, l_end, l_all, var_all, pl, -1);
    end
    %
    %% final disp
    %
    if Opt.display
        fprintf('Time Elapsed: %.3f s\n',toc(t_display));
    end
    %
end