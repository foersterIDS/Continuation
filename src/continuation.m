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
    [Opt,ds0,Opt_is_set] = continuation.input(varargin,fun,var0,l_start,l_end,ds0);
    [Opt,ds0,Stepsize_options] = step_size.initialize(Opt,var0,l_start,l_end,ds0);
    if Stepsize_options.rate_of_contraction
        global solver_stepsizes;
    end
    [Bifurcation,Counter,Do,Info,Initial,Path,Plot,Solver] = aux.initialize_structs(var0,l_start,l_end,ds0,Opt,Stepsize_options.rate_of_contraction);
    clear('var0','l_start','l_end','ds0');
    res_corr = continuation.corrector(Opt);
    ds = Info.ds0;
    aux.print_line(Opt,'Starting path continuation...\n');
    t_display = tic;
    %
    %% find initial solution
    %
    residual_initial = @(v) aux.residual_fixed_value(fun,v,Opt.l_0,Opt);
    [Path.var_all,fun_initial,initial_exitflag,solver_output,initial_jacobian] = Solver.main(residual_initial,Info.var0,Opt.dscale0(1:end-1));
    solver_jacobian = initial_jacobian;
    x_plus = [];
    if Stepsize_options.predictor
        x_predictor_plus = [];
    end
    last_jacobian = [];
    break_fun_out = [];
    if initial_exitflag>0
        Path.l_all = Opt.l_0;
        Path.s_all = 0;
        Do.continuation = true;
        Do.loop = true;
        solver_jacobian = [solver_jacobian,aux.numeric_jacobian(@(x) fun(x(1:Info.nv),x(Info.nv+1)),[Path.var_all;Opt.l_0],'central_value',fun_initial,'derivative_dimensions',Info.nv+1,'diffquot',Opt.diffquot)];
        previous_jacobian = solver_jacobian;
        last_jacobian = solver_jacobian;
        [~,break_fun_out] = Opt.break_function(fun_initial,solver_jacobian,Path.var_all,Path.l_all,break_fun_out);
        aux.print_line(Opt,'Initial solution at l = %.2e\n',Opt.l_0);
        if aux.ison(Opt.bifurcation)
            sign_det_jacobian = sign(det(initial_jacobian));
        end
        if aux.ison(Opt.plot)
            [Plot, Opt] = plot.live_plot(Opt, Info, Path, Info.ds0, Info.ds0, solver_output.iterations, Counter);
        end
    else
        Path.var_all = [];
        Path.l_all = [];
        Path.s_all = [];
        exitflag = -2; % no initial solution found
        Info.exitflag = -2;
        var_all = Path.var_all;
        l_all = Path.l_all;
        s_all = Path.s_all;
        Do.continuation = false;
        aux.print_line(Opt,'No initial solution found.\n');
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
        %% save old values of stepsize information
        %
        if Stepsize_options.iterations
            if isfield(solver_output, 'iterations')
                iterations_tmp = solver_output.iterations(end);
            else
                iterations_tmp = [];
            end
        end
        %
        if Stepsize_options.speed_of_continuation
            if isfield(Path, 'speed_of_continuation') && ~isempty(Path.speed_of_continuation)
                speed_of_continuation_tmp = Path.speed_of_continuation(end);
            else
                speed_of_continuation_tmp = [];
            end
        end
        %
        if Stepsize_options.predictor
            if isfield(Path, 'predictor') && ~isempty(Path.x_predictor)
                predictor_tmp = Path.x_predictor(:,end);
            else
                predictor_tmp = [];
            end
        end
        %
        if Stepsize_options.rate_of_contraction
            if isfield(solver_output, 'rate_of_contraction')
                rate_of_contraction_tmp = solver_output.rate_of_contraction;
            else
                rate_of_contraction_tmp = [];
            end
        end
        %
        %% residual
        %
        if Do.deflate
            try
                residual = @(x) aux.deflation(residual,x_deflation,x,Opt);
            catch
                aux.print_line(Opt,'---> delation: catch!\n');
                Counter.catch = Counter.catch + 1;
                if Counter.catch >= 3
                    aux.print_line(Opt,'--> Error in input! catch was used too often!\n');
                    break;
                end
            end
        elseif Do.suspend
            residual = @(v,l_fix) aux.residual_suspend_continuation(fun,v,l_fix,Opt);
        else
            residual = @(x) aux.merge_residuals(fun,res_corr,x,[Path.var_all;Path.l_all],ds,last_jacobian,Opt);
        end
        %
        %% predictor
        %
        if Stepsize_options.speed_of_continuation
            tic;
        end
        %
        try
            if Do.homotopy
                %% Homotopy
                x_last_step = [var_all(:,end);Path.l_all(end)];
                x_predictor = homotopy.h_continuation(residual,x_last_step,Opt);
                var_predictor = x_predictor(1:end-1);
                l_predictor = x_predictor(end);
            else
                [var_predictor,l_predictor,fun_predictor,s_predictor,ds] = continuation.predictor(Path,ds,last_jacobian,fun,res_corr,Solver,Opt);
                x_predictor = [var_predictor;l_predictor];
            end
        catch
            [var_predictor,l_predictor,fun_predictor,s_predictor,ds] = continuation.predictor(Path,ds,last_jacobian,fun,res_corr,Solver,Opt);
            x_predictor = [var_predictor;l_predictor];
            aux.print_line(Opt,'---> predictor: catch!\n');
            Counter.catch = Counter.catch + 1;
            if Counter.catch >= 3
                aux.print_line(Opt,'--> Error in input! catch was used too often!\n');
                break;
            end
        end
        %
        %% solve
        %
        try
            dscale = aux.get_dscale(Opt,Path);
            if sign(Path.l_all(end)-Opt.l_target)*sign(l_predictor-Opt.l_target)<=0
                %% try to converge to target
                residual_target = @(v) aux.residual_fixed_value(fun,v,Opt.l_target,Opt);
                var_predictor_ctt = (var_predictor - Path.var_all(:,end))*(abs(Opt.l_target-Path.l_all(end))/abs(l_predictor-Path.l_all(end)))+Path.var_all(:,end);
                [var_solution,fun_solution,solver_exitflag,solver_output,solver_jacobian] = Solver.main(residual_target,var_predictor_ctt,dscale(1:end-1));
                x_solution = [var_solution;Opt.l_target];
                Do.convergeToTarget = true;
            else
                %% regular solver
                %            
                if Opt.solver.fsolve && Opt.solver_force1it
                    [x_solution,fun_solution,solver_exitflag,solver_output,solver_jacobian] = Solver.main(residual,x_predictor,dscale);
                    if solver_output.iterations < 1
                        % perturbate initial solution by tolerance of
                        % solver
                        pert = Opt.solver_tol * ones(numel(x_predictor),1) / numel(x_predictor);
                        [x_solution,fun_solution,solver_exitflag,solver_output,solver_jacobian] = Solver.main(residual,x_predictor + pert,dscale);
                    end
                else
                    if Do.suspend
                        [v_solution,fun_solution,solver_exitflag,solver_output,solver_jacobian] = Solver.main(@(v) residual(v,x_predictor(end)),x_predictor(1:(end-1)),dscale(1:(end-1)));
                        x_solution = [v_solution;x_predictor(end)];
                    else
                        [x_solution,fun_solution,solver_exitflag,solver_output,solver_jacobian] = Solver.main(residual,x_predictor,dscale);
                    end
                end
                Do.convergeToTarget = false;
            end
            is_current_jacobian = true;
        catch
            x_solution = NaN(size(x_predictor));
            fun_solution = inf(size(x_predictor));
            solver_exitflag = -2;
            solver_output = Solver.default_output;
            Do.convergeToTarget = false;
            aux.print_line(Opt,'---> solve: catch!\n');
            Counter.catch = Counter.catch + 1;
            if Counter.catch >= 3
                aux.print_line(Opt,'--> Error in input! catch was used too often!\n');
                break;
            end
        end
        %
        %% adaptive corrector
        %
        [Do,Opt,corr_info] = corrector.adapt(Do,Opt,Path,solver_exitflag,solver_output,Solver,fun,x_predictor,dscale,last_jacobian,ds);
        %
        %% check result
        %
        [val,is_reverse,catch_flag,inv_poi_str,Do,Opt] = aux.validate_result(x_solution,x_plus,fun_solution,Path,ds,solver_output,solver_exitflag,solver_jacobian,last_jacobian,fun_predictor,s_predictor,Do,Bifurcation,Info,Counter,Plot,Opt);
        if val
            %% valid result
            if isempty(x_plus)
                Path.var_all = [Path.var_all,x_solution(1:end-1)];
                Path.l_all = [Path.l_all,x_solution(end)];
                Path.s_all = [Path.s_all,Path.s_all(end)+norm(x_solution-[Path.var_all(:,end-1);Path.l_all(end-1)])];
                if Stepsize_options.predictor
                    Path.x_predictor = [Path.x_predictor, x_predictor];
                end
                previous_jacobian = last_jacobian;
                last_jacobian = solver_jacobian;
                Counter.valid_stepback = 0;
            else
                Path.var_all = [Path.var_all,x_solution(1:end-1),x_plus(1:end-1)];
                Path.l_all = [Path.l_all,x_solution(end),x_plus(end)];
                Path.s_all = [Path.s_all,Path.s_all(end)+norm(x_solution-[Path.var_all(:,end-2);Path.l_all(end-2)])*[1,1]+norm(x_plus-x_solution)*[0,1]];
                if Stepsize_options.predictor
                    Path.x_predictor = [Path.x_predictor, x_predictor, x_predictor_plus];
                    x_predictor_plus = [];
                end
                x_plus = [];
                previous_jacobian = solver_jacobian;
                last_jacobian = plus_jacobian;
                plus_jacobian = [];
                Counter.valid_stepback = Counter.valid_stepback+1;
            end
            Do.deflate = false;
            Do.homotopy = false;
            Do.stepback = false;
            Do.suspend = false;
            Counter.error = 0;
            Counter.step = Counter.step + 1;
            if Do.convergeToTarget
                last_jacobian = aux.get_jacobian(fun,Path.var_all(:,end),Path.l_all(end),Opt);
            end
            if Do.remove
                if Path.s_all(end)>(s_rmv+ds_rmv)
                    Opt.ds_max = Initial.ds_max;
                    Do.remove = false;
                    Counter.remove = 0;
                end
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
                %% stepback
                if Counter.valid_stepback<Opt.stepback_error_counter
                    Do.stepback = true;
                    x_plus = [Path.var_all(:,end);Path.l_all(end)];
                    if Stepsize_options.predictor
                        x_predictor_plus = x_predictor;
                    end
                else
                    Do.stepback = false;
                end
                Path.var_all(:,end) = [];
                Path.l_all(end) = [];
                Path.s_all(end) = [];
                if Stepsize_options.predictor
                    Path.x_predictor(:,end) = [];
                end
                plus_jacobian = last_jacobian;
                last_jacobian = previous_jacobian;
            elseif (Counter.error==Opt.stepback_error_counter+1) && (numel(Path.l_all)>1)
                %% undo stepback
                if ~isempty(x_plus)
                    Path.var_all = [Path.var_all,x_plus(1:end-1)];
                    Path.l_all = [Path.l_all,x_plus(end)];
                    Path.s_all = [Path.s_all,Path.s_all(end)+norm([Path.var_all(:,end);Path.l_all(end)]-[Path.var_all(:,end-1);Path.l_all(end-1)])];
                    if Stepsize_options.predictor
                        Path.x_predictor = [Path.x_predictor, x_predictor_plus];
                        x_predictor_plus = [];
                    end
                    x_plus = [];
                end
                Do.stepback = false;
                Do.suspend = false;
            elseif (Counter.error==Opt.suspend_continuation_error_counter) && (numel(Path.l_all)>1)
                %% suspend
                x_plus = [];
                if Stepsize_options.predictor
                    x_predictor_plus = [];
                end
                Do.stepback = false;
                Do.suspend = true;
            elseif logical(Opt.remove_error_counter) && ((Counter.error==Opt.remove_error_counter) && (numel(Path.l_all)>1))
                %% remove
                n_path = numel(Path.l_all);
                s_rmv = Path.s_all(n_path);
                n_rmv = min([2*Opt.remove_error_counter,n_path-1]);
                Opt.ds_max = max([mean(diff(Path.s_all(n_path+((-ceil(n_rmv/2)+1):0))))*0.75,Opt.ds_min]);
                Path.var_all(:,n_path+((-n_rmv+1):0)) = [];
                Path.l_all(n_path+((-n_rmv+1):0)) = [];
                Path.s_all(n_path+((-n_rmv+1):0)) = [];
                ds_rmv = s_rmv-Path.s_all(end);
                dscale_rmv = aux.get_dscale(Opt,Path);
                residual_fixed_value_rmv = @(v) aux.residual_fixed_value(fun,v,Path.l_all(end),Opt);
                [var_rmv,fun_rmv,~,~,rmv_jacobian] = Solver.main(residual_fixed_value_rmv,Path.var_all(:,end),dscale_rmv(1:end-1));
                if ~isempty(Bifurcation.bif) && numel(Bifurcation.bif(1,:))>0
                    n_bifs_rmv = sum(sum(Bifurcation.bif(1,:)'==(n_path+((-n_rmv+1):0))));
                    if n_bifs_rmv>0
                        Bifurcation.bif(:,end+((1-n_bifs_rmv):0)) = [];
                        sign_det_jacobian = sign_det_jacobian*(-1)^(n_bifs_rmv);
                    end
                end
                Path.var_all(:,end) = var_rmv;
                last_jacobian = [rmv_jacobian,aux.numeric_jacobian(@(x) fun(x(1:Info.nv),x(Info.nv+1)),[Path.var_all(:,end);Path.l_all(end)],'central_value',fun_rmv,'derivative_dimensions',Info.nv+1,'diffquot',Opt.diffquot)];
                x_plus = [];
                if Stepsize_options.predictor
                    x_predictor_plus = [];
                end
                Do.stepback = false;
                Do.suspend = false;
                Do.remove = true;
                Counter.remove = Counter.remove+1;
            else
                %% else
                x_plus = [];
                if Stepsize_options.predictor
                    x_predictor_plus = [];
                end
                Do.stepback = false;
                Do.suspend = false;
                Do.remove = false;
            end
            if Opt.include_reverse && is_reverse && solver_exitflag>0
                Path = aux.include_reverse(x_solution,Path);
            end
            if aux.ison(Opt.homotopy) && Counter.error>=Opt.homotopy_error_counter
                Do.homotopy = true;
            else
                Do.homotopy = false;
            end
            if catch_flag
                Counter.catch = Counter.catch + 1;
                if Counter.catch >= 3
                    aux.print_line(Opt,'--> Error in input! catch was used too often!\n');
                    break;
                end
            end
        end
        %
        %% Bifurcations
        %
        Bifurcation.flag = 0;
        if aux.ison(Opt.bifurcation) && val && ~Do.homotopy && numel(Path.l_all)>2
            if ~is_current_jacobian
                %% get jacobian if not current
                solver_jacobian = aux.get_jacobian(fun,Path.var_all(:,end),Path.l_all(end),Opt);
            end
            [Bifurcation,sign_det_jacobian,Path] = bifurcation.check(fun,solver_jacobian(1:Info.nv,1:Info.nv),Path,Bifurcation,sign_det_jacobian,res_corr,Solver,Opt);
        elseif aux.ison(Opt.bifurcation) && val && numel(Path.l_all)<=2
            if ~is_current_jacobian
                %% get jacobian if not current
                solver_jacobian = aux.get_jacobian(fun,Path.var_all(:,end),Path.l_all(end),Opt);
            end
            sign_det_jacobian = sign(det(solver_jacobian(1:Info.nv,1:Info.nv)));
        end
        %
        %% step size control
        %
        % save latest stepsize
        %
        dsim1 = ds;
        %
        % update stepsize information
        % measure speed
        %
        if Stepsize_options.speed_of_continuation
            time_needed = toc;
            if ~isempty(speed_of_continuation_tmp)
                Path.speed_of_continuation = [Path.speed_of_continuation(end), ds/time_needed];
            else
                Path.speed_of_continuation = ds/time_needed;
            end
            
        end
        %
        % measure rate of contraction
        if Stepsize_options.rate_of_contraction
            if size(solver_stepsizes, 1) < 3
                if ~isempty(rate_of_contraction_tmp)
                    solver_output.rate_of_contraction = [rate_of_contraction_tmp(end), Opt.optimal_contraction_rate];
                else
                    solver_output.rate_of_contraction = Opt.optimal_contraction_rate;
                end
            else
                if ~isempty(rate_of_contraction_tmp)
                    solver_output.rate_of_contraction = [rate_of_contraction_tmp(end), solver_stepsizes(3,2)/solver_stepsizes(2,2)];
                else
                    solver_output.rate_of_contraction = solver_stepsizes(3,2)/solver_stepsizes(2,2);
                end
            end
        end
        %
        if Stepsize_options.iterations
            if ~isempty(iterations_tmp)
                solver_output.iterations = [iterations_tmp, solver_output.iterations];
            end
        end
        %
        if Stepsize_options.predictor
            if ~isempty(predictor_tmp)
                Path.predictor = [predictor_tmp, Path.x_predictor];
            end
        end
        %
        % adjust stepsize
        %
        ds = step_size.control(ds,Counter,solver_output,Do,x_plus,Path,last_jacobian,Opt);
        %
        %% end loop
        %
        if val
            aux.print_line(Opt,'-----> continued at l = %.4e\t|\tnew step size: ds = %.2e\t|\tloop counter = %d\t|\tstep = %d\t|\titerations = %d/%d\n',Path.l_all(end),ds,Counter.loop,Counter.step,solver_output.iterations(end),Opt.n_iter_opt);
        else
            aux.print_line(Opt,'-----> invalid point %s |\tnew step size: ds = %.2e\t|\tloop counter = %d\t|\tstep = %d\t|\titerations = %d/%d\n',inv_poi_str,ds,Counter.loop,Counter.step,solver_output.iterations(end),Opt.n_iter_opt);
        end
        [Do,Info,Path,break_fun_out,Opt] = aux.exit_loop(Do,Info,Path,Opt,Counter,Bifurcation,ds,fun_solution,solver_jacobian,break_fun_out,val);
        exitflag = Info.exitflag;
        var_all = Path.var_all;
        l_all = Path.l_all;
        s_all = Path.s_all;
        if Do.change_corrector
            Opt = aux.seton(Opt,'corrector',corr_info);
            res_corr = continuation.corrector(Opt);
            Do.change_corrector = false;
        end
        Opt = aux.update_Opt(Opt,Opt_is_set,Info);
        %
        %% live plot
        %
        if aux.ison(Opt.plot) && val
            try
                [Plot, Opt] = plot.live_plot(Opt, Info, Path, ds, dsim1, solver_output.iterations, Counter, fun_predictor, s_predictor, Plot, Bifurcation);
            catch
                aux.print_line(Opt,'--> The plot update has failed.\n');
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
            delete(Plot.pl_curr);
            [Path,Bifurcation] = bifurcation.trace(Opt,Path,Bifurcation,Solver,Info,fun,res_corr);
            last_jacobian = [];
        catch
            aux.print_line(Opt,'--> Failed to trace bifurcations.\n');
        end
    end
    %
    %% live plot finalization
    %
    if aux.ison(Opt.plot) && initial_exitflag>0
        try
            Bifurcation_last_plot = Bifurcation;
            Bifurcation_last_plot.flag = -1;
            plot.live_plot(Opt, Info, Path, ds, dsim1, solver_output.iterations, Counter, fun_predictor, s_predictor, Plot, Bifurcation_last_plot);
            if isfield(Plot,'pl_curr')
                delete(Plot.pl_curr);
            end
        catch
            aux.print_line(Opt,'--> The plot update has failed.\n');
        end
    end
    %
    %% final disp
    %
    aux.print_line(Opt,[Info.exit_msg,'\n']);
    aux.print_line(Opt,'--> time elapsed: %.3f s\n',toc(t_display));
    %
end