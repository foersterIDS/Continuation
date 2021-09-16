%% path continuation - continuation
%
%   Institute of Dynamics and Vibration Research
%   Leibniz University Hannover
%   08.05.2020 - Alwin Förster
%
%   [var_all,l_all,exitflag,Bifurcation] = continuation(fun,var0,l_start,l_end,ds0,varargin)
%
%   fun = fun(var,l) != 0
%   l_start <= l <= l_end
%   ds0: initial stepsize
%
function [var_all,l_all,exitflag,Bifurcation,s_all,last_jacobian,break_fun_out] = continuation(fun,var0,l_start,l_end,ds0,varargin)
    %% initialize
    %
    warning on;
    [Opt,ds0] = continuation_input(varargin,fun,var0,l_start,l_end,ds0);
    [Bifurcation,Counter,Do,Info,Path,Plot] = initialize_structs(var0,l_start,l_end);
    [solver,predictor_solver,default_solver_output] = continuation_solver(Opt);
    res_corr = residual_corrector(Opt);
    ds = ds0;
    if Opt.display
        fprintf('Starting path continuation...\n');
        t_display = tic;
    end
    %
    %% find initial solution
    %
    residual_initial = @(v) residual_fixed_value(fun,v,Opt.l_0,Opt);
    [Path.var_all,fun_initial,initial_exitflag,solver_output,initial_jacobian] = solver(residual_initial,var0,Opt.dscale0(1:end-1));
    solver_jacobian = initial_jacobian;
    x_plus = [];
    last_jacobian = [];
    break_fun_out = [];
    if initial_exitflag>0
        Path.l_all = Opt.l_0;
        Path.s_all = 0;
        Do.continuation = true;
        Do.loop = true;
        solver_jacobian = [solver_jacobian,numeric_jacobian(@(x) fun(x(1:Info.nv),x(Info.nv+1)),[Path.var_all;Opt.l_0],'central_value',fun_initial,'derivative_dimensions',Info.nv+1,'diffquot',Opt.diffquot)];
        previous_jacobian = solver_jacobian;
        last_jacobian = solver_jacobian;
        [~,break_fun_out] = Opt.break_function(fun_initial,solver_jacobian,Path.var_all,Path.l_all,break_fun_out);
        if Opt.display
            fprintf('Initial solution at l = %.2e\n',Opt.l_0);
        end
        if ison(Opt.bifurcation)
            sign_det_jacobian = sign(det(initial_jacobian));
        end
        if ison(Opt.plot)
            [Plot, Opt] = live_plot(Opt, Info, Path, ds0, ds0, solver_output.iterations, Counter);
        end
    else
        Path.var_all = [];
        Path.l_all = [];
        Path.s_all = [];
        exitflag = -2; % no initial solution found
        Info.exitflag = -2;
        Do.continuation = false;
        if Opt.display
            fprintf('No initial solution found.\n');
        end
    end
    %
    %% continuation
    %
    while Do.continuation
        %% initialize loop
        %
        Counter.loop = Counter.loop+1;
        is_current_jacobian = false;
        Counter.catch_old = Counter.catch;
        %
        %% residual
        %
        residual = @(x) merge_residuals(fun,res_corr,x,[Path.var_all;Path.l_all],ds,Opt);
        if Do.deflate
            try
                residual = @(x) deflation(residual,x_deflation,x,Opt);
            catch
                error('Error occured during deflation.');
                Counter.catch = Counter.catch + 1;
                if Counter.catch >= 3
                    warning('Error in input! catch was used too often!');
                    break;
                end
            end
        end
        %
        %% predictor
        %
        try
            if Do.homotopy
                %% Homotopy
                x_last_step = [var_all(:,end);Path.l_all(end)];
                x_predictor = homotopy(residual,x_last_step,Opt);
                var_predictor = x_predictor(1:end-1);
                l_predictor = x_predictor(end);
            else
                [var_predictor,l_predictor,fun_predictor,s_predictor] = predictor(Path,ds,last_jacobian,fun,res_corr,predictor_solver,Opt);
                x_predictor = [var_predictor;l_predictor];
            end
        catch
            [var_predictor,l_predictor,fun_predictor,s_predictor] = predictor(Path,ds,last_jacobian,fun,res_corr,predictor_solver,Opt);
            x_predictor = [var_predictor;l_predictor];
            warning('predictor: catch!');
            Counter.catch = Counter.catch + 1;
            if Counter.catch >= 3
                    warning('Error in input! catch was used too often!');
                    break;
            end
        end
        %
        %% solve
        %
        try
            dscale = get_dscale(Opt,Path);
            if sign(Path.l_all(end)-Opt.l_target)*sign(l_predictor-Opt.l_target)<=0
                %% try to converge to target
                residual_target = @(v) residual_fixed_value(fun,v,Opt.l_target,Opt);
                var_predictor_ctt = (var_predictor - Path.var_all(:,end))*(abs(Opt.l_target-Path.l_all(end))/abs(l_predictor-Path.l_all(end)))+Path.var_all(:,end);
                [var_solution,fun_solution,solver_exitflag,solver_output,solver_jacobian] = solver(residual_target,var_predictor_ctt,dscale(1:end-1));
                x_solution = [var_solution;Opt.l_target];
                Do.convergeToTarget = true;
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
                Do.convergeToTarget = false;
            end
            is_current_jacobian = true;
        catch
            x_solution = NaN(size(x_predictor));
            fun_solution = inf(size(x_predictor));
            solver_exitflag = -2;
            solver_output = default_solver_output;
            Do.convergeToTarget = false;
            warning('solve: catch!');
            Counter.catch = Counter.catch + 1;
            if Counter.catch >= 3
                    warning('Error in input! catch was used too often!');
                    break;
            end
        end
        %
        %% check result
        %
        [val,is_reverse,catch_flag,Do,Opt] = validate_result(x_solution,x_plus,fun_solution,Path,ds,solver_output,solver_exitflag,solver_jacobian,last_jacobian,fun_predictor,s_predictor,Do,Bifurcation,Info,Counter,Plot,Opt);
        if val
            %% valid result
            if isempty(x_plus)
                Path.var_all = [Path.var_all,x_solution(1:end-1)];
                Path.l_all = [Path.l_all,x_solution(end)];
                Path.s_all = [Path.s_all,Path.s_all(end)+norm(x_solution-[Path.var_all(:,end-1);Path.l_all(end-1)])];
                previous_jacobian = last_jacobian;
                last_jacobian = solver_jacobian;
                Counter.valid_stepback = 0;
            else
                Path.var_all = [Path.var_all,x_solution(1:end-1),x_plus(1:end-1)];
                Path.l_all = [Path.l_all,x_solution(end),x_plus(end)];
                Path.s_all = [Path.s_all,Path.s_all(end)+norm(x_solution-[Path.var_all(:,end-2);Path.l_all(end-2)])*[1,1]+norm(x_plus-x_solution)*[0,1]];
                x_plus = [];
                previous_jacobian = solver_jacobian;
                last_jacobian = plus_jacobian;
                plus_jacobian = [];
                Counter.valid_stepback = Counter.valid_stepback+1;
            end
            Do.deflate = false;
            Do.homotopy = false;
            Do.stepback = false;
            Counter.error = 0;
            Counter.step = Counter.step + 1;
            if Do.convergeToTarget
                last_jacobian = get_jacobian(fun,Path.var_all(:,end),Path.l_all(end),Opt);
            end
        else
            %% invalid result
            Counter.error = Counter.error+1;
            if Opt.deflation && ~Do.deflate && ~isnan(sum(x_solution(:,end))) && is_reverse && Counter.error>=Opt.deflation_error_counter
                Do.deflate = true;
                x_deflation = x_solution;
            else
                Do.deflate = false;
            end
            if (Counter.error==Opt.stepback_error_counter) && (numel(Path.l_all)>1)
                if Counter.valid_stepback<Opt.stepback_error_counter
                    Do.stepback = true;
                    x_plus = [Path.var_all(:,end);Path.l_all(end)];
                else
                    Do.stepback = false;
                end
                Path.var_all(:,end) = [];
                Path.l_all(end) = [];
                Path.s_all(end) = [];
                plus_jacobian = last_jacobian;
                last_jacobian = previous_jacobian;
            elseif (Counter.error==Opt.stepback_error_counter+1) && (numel(Path.l_all)>1)
                Path.var_all = [Path.var_all,x_plus(1:end-1)];
                Path.l_all = [Path.l_all,x_plus(end)];
                Path.s_all = [Path.s_all,Path.s_all(end)+norm([Path.var_all(:,end);Path.l_all(end)]-[Path.var_all(:,end-1);Path.l_all(end-1)])];
                x_plus = [];
                Do.stepback = false;
            else
                x_plus = [];
                Do.stepback = false;
            end
            if Opt.include_reverse && is_reverse && solver_exitflag>0
                Path = include_reverse(x_solution,Path);
            end
            if ison(Opt.homotopy) && Counter.error>=Opt.homotopy_error_counter
                Do.homotopy = true;
            else
                Do.homotopy = false;
            end
            if catch_flag
                Counter.catch = Counter.catch + 1;
                if Counter.catch >= 3
                    warning('Error in input! catch was used too often!');
                    break;
                end
            end
        end
        %
        %% Bifurcations
        %
        Bifurcation.flag = 0;
        if ison(Opt.bifurcation) && val && ~Do.homotopy && numel(Path.l_all)>2
            if ~is_current_jacobian
                %% get jacobian if not current
                solver_jacobian = get_jacobian(fun,Path.var_all(:,end),Path.l_all(end),Opt);
            end
            [Bifurcation,sign_det_jacobian,Path] = check_bifurcation(fun,solver_jacobian(1:Info.nv,1:Info.nv),Path,Bifurcation,sign_det_jacobian,res_corr,predictor_solver,Opt);
        elseif ison(Opt.bifurcation) && val && numel(Path.l_all)<=2
            if ~is_current_jacobian
                %% get jacobian if not current
                solver_jacobian = get_jacobian(fun,Path.var_all(:,end),Path.l_all(end),Opt);
            end
            sign_det_jacobian = sign(det(solver_jacobian(1:Info.nv,1:Info.nv)));
        end
        %
        %% step size control
        %
        dsim1 = ds;
        ds = step_size_control(ds,ds0,Counter,solver_output,Do,x_plus,Path,Opt);
        %
        %% end loop
        %
        if Opt.display
            if val
                fprintf('-----> continued at l = %.4e\t|\tnew arc-length: ds = %.2e\t|\tloop counter = %d\t|\tstep = %d\t|\titerations = %d/%d\n',Path.l_all(end),ds,Counter.loop,Counter.step,solver_output.iterations,Opt.n_iter_opt);
            else
                fprintf('-----> invalid point\t\t\t\t|\tnew arc-length: ds = %.2e\t|\tloop counter = %d\t|\tstep = %d\t|\titerations = %d/%d\n',ds,Counter.loop,Counter.step,solver_output.iterations,Opt.n_iter_opt);
            end
        end
        [Do,Info,Path,break_fun_out,Opt] = exit_loop(Do,Info,l_start,l_end,Path,Opt,Counter,Bifurcation,ds,fun_solution,solver_jacobian,break_fun_out,val);
        exitflag = Info.exitflag;
        var_all = Path.var_all;
        l_all = Path.l_all;
        s_all = Path.s_all;
        %
        %% live plot
        %
        if ison(Opt.plot) && val
            try
                [Plot, Opt] = live_plot(Opt, Info, Path, ds, dsim1, solver_output.iterations, Counter, fun_predictor, s_predictor, Plot, Bifurcation);
            catch
                warning('The plot update has failed.');
            end
        end
        %
        %% catch counter
        %
        if Counter.catch_old == Counter.catch
            Counter.catch = 0;
        end
        %
    end
    %
    %% bifurcation tracing
    %
    if Opt.bifurcation.trace
        try
            [Path,Bifurcation] = trace_bifurcations(Opt,Path,Bifurcation,solver,fun,l_start,l_end,res_corr,predictor_solver);
            last_jacobian = [];
        catch
            warning('Failed to trace bifurcations.');
        end
    end
    %
    %% live plot finalization
    %
    if ison(Opt.plot) && initial_exitflag>0
        try
            Bifurcation_last_plot = Bifurcation;
            Bifurcation_last_plot.flag = -1;
            live_plot(Opt, Info, Path, ds, dsim1, solver_output.iterations, Counter, fun_predictor, s_predictor, Plot, Bifurcation_last_plot);
        catch
            warning('The plot update has failed.');
        end
    end
    %
    %% final disp
    %
    if Opt.display
        fprintf('Time Elapsed: %.3f s\n',toc(t_display));
    end
    %
end