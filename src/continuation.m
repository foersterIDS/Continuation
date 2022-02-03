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
function [var_all,l_all,exitflag,Bifurcation,s_all,jacobian_out,break_fun_out,Info_out] = ...
    continuation(fun,var0,l_start,l_end,ds0,varargin)
    %% initialize
    %
    warning on;
    [Opt,ds0,Opt_is_set] = continuation.input(varargin,fun,var0,l_start,l_end,ds0);
    [Opt,ds0,Stepsize_options] = step_size.initialize(Opt,var0,l_start,l_end,ds0);
    if Stepsize_options.rate_of_contraction
        global solver_stepsizes;
    end
    [Bifurcation,Counter,Do,Info,Info_out,Initial,Jacobian,Path,Plot,Plus,Solver,Temp] = aux.initialize_structs(var0,l_start,l_end,ds0,Opt,Stepsize_options.rate_of_contraction);
    clear('var0','l_start','l_end','ds0');
    res_corr = continuation.corrector(Opt);
    ds = Info.ds0;
    aux.print_line(Opt,'Starting path continuation...\n');
    t_display = tic;
    %
    %% find initial solution
    %
    residual_initial = @(v) aux.residual_fixed_value(fun,v,Opt.l_0,Opt);
    [Path.var_all,fun_initial,initial_exitflag,Solver.output,Jacobian.initial] = Solver.main(residual_initial,Info.var0,Opt.dscale0(1:end-1));
    Jacobian.solver = Jacobian.initial;
    break_fun_out = [];
    event_out = false;
    if initial_exitflag>0
        Path.l_all = Opt.l_0;
        Path.s_all = 0;
        Do.continuation = true;
        Do.loop = true;
        Jacobian.solver = [Jacobian.solver,aux.numeric_jacobian(@(x) fun(x(1:Info.nv),x(Info.nv+1)),[Path.var_all;Opt.l_0],'central_value',fun_initial,'derivative_dimensions',Info.nv+1,'diffquot',Opt.diffquot)];
        Jacobian.previous = Jacobian.solver;
        Jacobian.last = Jacobian.solver;
        [~,break_fun_out] = Opt.break_function(fun_initial,Jacobian.solver,Path.var_all,Path.l_all,break_fun_out);
        if Opt.step_size_event
            [ds,Counter,event_out,~] = step_size.event_adjustment(ds,Path,Counter,Opt,event_out);
        end
        aux.print_line(Opt,'Initial solution at l = %.2e\n',Opt.l_0);
        if aux.ison(Opt.bifurcation)
            Jacobian.sign_det = sign(det(Jacobian.initial));
        end
        if aux.ison(Opt.plot)
            [Plot, Opt] = plot.live_plot(Opt, Info, Path, Info.ds0, Info.ds0, Solver.output.iterations(end), Counter);
        end
        if Stepsize_options.rate_of_contraction
            if size(solver_stepsizes, 1) < 3
                Solver.output.rate_of_contraction = Opt.optimal_contraction_rate;
            else
                Solver.output.rate_of_contraction = solver_stepsizes(3,2)/solver_stepsizes(2,2);
            end
        end
        if Stepsize_options.speed_of_continuation
             Path.speed_of_continuation = Opt.speed_of_continuation;
        end
        %
        if Stepsize_options.predictor
            Path.x_predictor = [Info.var0;Info.l_start];
        end
    else
        Path.var_all = [];
        Path.l_all = [];
        Path.s_all = [];
        Info.exitflag = -2;
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
            Temp.length = length(Solver.output.iterations);
            Temp.flag = 0;
            if Temp.length < 3 && Temp.length > 0
                to_take = 1:Temp.length;
            elseif Temp.length == 0
                Temp.flag = 1;
            else
                to_take = Temp.length-1:Temp.length;
            end
            if isfield(Solver.output, 'iterations') && ~isempty(Solver.output.iterations) && ~Temp.flag
                Temp.iterations = Solver.output.iterations(to_take);
            else
                Temp.iterations = [];
            end
        end
        %
        if Stepsize_options.speed_of_continuation
            Temp.length = length(Path.speed_of_continuation);
            Temp.flag = 0;
            if Temp.length < 3 && Temp.length > 0
                to_take = 1:Temp.length;
            elseif Temp.length == 0
                Temp.flag = 1;
            else
                to_take = Temp.length-1:Temp.length;
            end
            if isfield(Path, 'speed_of_continuation') && ~isempty(Path.speed_of_continuation) && ~Temp.flag
                Temp.speed_of_continuation = Path.speed_of_continuation(to_take);
            else
                Temp.speed_of_continuation = [];
            end
        end
        %
        if Stepsize_options.predictor
            Temp.length = size(Path.x_predictor,2);
            Temp.flag = 0;
            if Temp.length < 3 && Temp.length > 0
                to_take = 1:Temp.length;
            elseif Temp.length == 0
                Temp.flag = 1;
            else
                to_take = Temp.length-1:Temp.length;
            end
            if isfield(Path, 'x_predictor') && ~isempty(Path.x_predictor) && ~Temp.flag
                Temp.predictor = Path.x_predictor(:,to_take);
            else
                Temp.predictor = [];
            end
        end
        %
        if Stepsize_options.rate_of_contraction
            Temp.length = length(Solver.output.rate_of_contraction);
            Temp.flag = 0;
            if Temp.length < 3 && Temp.length > 0
                to_take = 1:Temp.length;
            elseif Temp.length == 0
                Temp.flag = 1;
            else
                to_take = Temp.length-1:Temp.length;
            end
            if isfield(Solver.output, 'rate_of_contraction') && ~isempty(Solver.output.rate_of_contraction) && ~Temp.flag
                Temp.rate_of_contraction = Solver.output.rate_of_contraction(to_take);
            else
                Temp.rate_of_contraction = [];
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
            residual = @(x) aux.merge_residuals(fun,res_corr,x,[Path.var_all;Path.l_all],ds,Jacobian.last,Opt);
        end
        %
        %% predictor
        %
        if Stepsize_options.speed_of_continuation
            time_needed = tic;
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
                [var_predictor,l_predictor,fun_predictor,s_predictor,ds] = continuation.predictor(Path,ds,Jacobian.last,fun,res_corr,Solver,Opt);
                x_predictor = [var_predictor;l_predictor];
            end
        catch
            [var_predictor,l_predictor,fun_predictor,s_predictor,ds] = continuation.predictor(Path,ds,Jacobian.last,fun,res_corr,Solver,Opt);
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
                [var_solution,fun_solution,Solver.exitflag,Solver.output,Jacobian.solver] = Solver.main(residual_target,var_predictor_ctt,dscale(1:end-1));
                x_solution = [var_solution;Opt.l_target];
                Do.convergeToTarget = true;
            else
                %% regular solver
                %            
                if Opt.solver.fsolve && Opt.solver_force1it
                    [x_solution,fun_solution,Solver.exitflag,Solver.output,Jacobian.solver] = Solver.main(residual,x_predictor,dscale);
                    if Solver.output.iterations(end) < 1
                        % perturbate initial solution by tolerance of
                        % solver
                        pert = Opt.solver_tol * ones(numel(x_predictor),1) / numel(x_predictor);
                        [x_solution,fun_solution,Solver.exitflag,Solver.output,Jacobian.solver] = Solver.main(residual,x_predictor + pert,dscale);
                    end
                else
                    if Do.suspend
                        [v_solution,fun_solution,Solver.exitflag,Solver.output,Jacobian.solver] = Solver.main(@(v) residual(v,x_predictor(end)),x_predictor(1:(end-1)),dscale(1:(end-1)));
                        x_solution = [v_solution;x_predictor(end)];
                    else
                        [x_solution,fun_solution,Solver.exitflag,Solver.output,Jacobian.solver] = Solver.main(residual,x_predictor,dscale);
                    end
                end
                Do.convergeToTarget = false;
            end
            is_current_jacobian = true;
        catch
            x_solution = NaN(size(x_predictor));
            fun_solution = inf(size(x_predictor));
            Solver.exitflag = -2;
            Solver.output = Solver.default_output;
            Do.convergeToTarget = false;
            aux.print_line(Opt,'---> solve: catch!\n');
            Counter.catch = Counter.catch + 1;
            if Counter.catch >= 3
                aux.print_line(Opt,'--> Error in input! catch was used too often!\n');
                break;
            end
        end
        %
        % calc stepsize information
        % measure speed
        %
        if Stepsize_options.speed_of_continuation
            time_needed = toc(time_needed);
            speed_of_continuation = ds/time_needed; 
        end
        %
        % measure rate of contraction
        if Stepsize_options.rate_of_contraction
            if size(solver_stepsizes, 1) < 3
                rate_of_contraction = Opt.optimal_contraction_rate;
            else
                rate_of_contraction = solver_stepsizes(3,2)/solver_stepsizes(2,2);
            end
        end
        %
        if Stepsize_options.iterations
            if ~isempty(Temp.iterations)
                Temp.iterations = [Temp.iterations, Solver.output.iterations];
            else
                Temp.iterations = Solver.output.iterations;
            end
        end
        %% adaptive corrector
        %
        [Do,Opt,corr_info] = corrector.adapt(Do,Opt,Path,Solver,fun,x_predictor,dscale,Jacobian.last,ds);
        %
        %% check result
        %
        [val,is_reverse,catch_flag,inv_poi_str,Do,Opt] = aux.validate_result(x_solution,Plus,fun_solution,Path,ds,Solver,Jacobian,fun_predictor,s_predictor,Do,Bifurcation,Info,Counter,Plot,Opt);
        if val
            %% valid result
            if isempty(Plus.x)
                Path.var_all = [Path.var_all,x_solution(1:end-1)];
                Path.l_all = [Path.l_all,x_solution(end)];
                Path.s_all = [Path.s_all,Path.s_all(end)+norm(x_solution-[Path.var_all(:,end-1);Path.l_all(end-1)])];
                % update stepsize information
                % measure speed
                %
                if Stepsize_options.speed_of_continuation
                    if ~isempty(Temp.speed_of_continuation)
                        Temp.speed_of_continuation = [Temp.speed_of_continuation, speed_of_continuation];
                    else
                        Temp.speed_of_continuation = speed_of_continuation;
                    end
                end
                %
                % measure rate of contraction
                if Stepsize_options.rate_of_contraction
                    if ~isempty(Temp.rate_of_contraction)
                        Temp.rate_of_contraction = [Temp.rate_of_contraction, rate_of_contraction];
                    else
                        Temp.rate_of_contraction = rate_of_contraction;
                    end
                end
                %
                if Stepsize_options.predictor
                    if ~isempty(Temp.predictor)
                        Temp.predictor = [Temp.predictor, x_predictor];
                    end
                end
                Jacobian.previous = Jacobian.last;
                Jacobian.last = Jacobian.solver;
                Counter.valid_stepback = 0;
            else
                Path.var_all = [Path.var_all,x_solution(1:end-1),Plus.x(1:end-1)];
                Path.l_all = [Path.l_all,x_solution(end),Plus.x(end)];
                Path.s_all = [Path.s_all,Path.s_all(end)+norm(x_solution-[Path.var_all(:,end-2);Path.l_all(end-2)])*[1,1]+norm(Plus.x-x_solution)*[0,1]];
                if Stepsize_options.predictor
                    Temp.predictor = [Temp.predictor, x_predictor, Plus.x_predictor];
                    Plus.x_predictor = [];
                end
                if Stepsize_options.speed_of_continuation
                    Temp.speed_of_continuation = [Temp.speed_of_continuation, speed_of_continuation, Plus.speed_of_continuation];
                    Plus.speed_of_continuation = [];
                end
                if Stepsize_options.rate_of_contraction
                    Temp.rate_of_contraction = [Temp.rate_of_contraction, rate_of_contraction, Plus.rate_of_contraction];
                    Plus.rate_of_contraction = [];
                end
                Plus.x = [];
                Jacobian.previous = Jacobian.solver;
                Jacobian.last = Plus.jacobian;
                Plus.jacobian = [];
                Counter.valid_stepback = Counter.valid_stepback+1;
            end
            Do.deflate = false;
            Do.homotopy = false;
            Do.stepback = false;
            Do.suspend = false;
            Counter.error = 0;
            Counter.step = Counter.step + 1;
            if Do.convergeToTarget
                Jacobian.last = aux.get_jacobian(fun,Path.var_all(:,end),Path.l_all(end),Opt);
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
                    Plus.x = [Path.var_all(:,end);Path.l_all(end)];
                    if Stepsize_options.predictor
                        Plus.x_predictor = x_predictor;
                    end
                    if Stepsize_options.speed_of_continuation
                        Plus.speed_of_continuation = speed_of_continuation;
                    end
                    if Stepsize_options.rate_of_contraction
                        Plus.rate_of_contraction = rate_of_contraction;
                    end
                else
                    Do.stepback = false;
                end
                Path.var_all(:,end) = [];
                Path.l_all(end) = [];
                Path.s_all(end) = [];
                if Stepsize_options.predictor 
                    if ~isempty(Temp.predictor)
                        Temp.predictor(:,end) = [];
                    else
                        Temp.predictor = [];
                    end
                end
                if Stepsize_options.speed_of_continuation
                    Temp.speed_of_continuation(end) = [];
                end
                if Stepsize_options.rate_of_contraction
                    Temp.rate_of_contraction(end) = [];
                end
                Plus.jacobian = Jacobian.last;
                Jacobian.last = Jacobian.previous;
            elseif (Counter.error==Opt.stepback_error_counter+1) && (numel(Path.l_all)>1)
                %% undo stepback
                if ~isempty(Plus.x)
                    Path.var_all = [Path.var_all,Plus.x(1:end-1)];
                    Path.l_all = [Path.l_all,Plus.x(end)];
                    Path.s_all = [Path.s_all,Path.s_all(end)+norm([Path.var_all(:,end);Path.l_all(end)]-[Path.var_all(:,end-1);Path.l_all(end-1)])];
                    if Stepsize_options.predictor
                        Temp.predictor = [Temp.predictor, Plus.x_predictor];
                    end
                    if Stepsize_options.speed_of_continuation
                        Temp.speed_of_continuation = [Temp.speed_of_continuation, Plus.speed_of_continuation];
                    end
                    if Stepsize_options.rate_of_contraction
                        Temp.rate_of_contraction = [Temp.rate_of_contraction, Plus.rate_of_contraction];
                    end
                    Plus = aux.clear_struct(Plus);
                end
                Do.stepback = false;
                Do.suspend = false;
            elseif (Counter.error==Opt.suspend_continuation_error_counter) && (numel(Path.l_all)>1)
                %% suspend
                Plus = aux.clear_struct(Plus);
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
                        Jacobian.sign_det = Jacobian.sign_det*(-1)^(n_bifs_rmv);
                    end
                end
                Path.var_all(:,end) = var_rmv;
                Jacobian.last = [rmv_jacobian,aux.numeric_jacobian(@(x) fun(x(1:Info.nv),x(Info.nv+1)),[Path.var_all(:,end);Path.l_all(end)],'central_value',fun_rmv,'derivative_dimensions',Info.nv+1,'diffquot',Opt.diffquot)];
                if Stepsize_options.predictor
                    if n_rmv >= 3
                        Temp.predictor = x_predictor;
                    else
                        Temp.predictor = Temp.predictor(:,1:end-n_rmv);
                    end
                end
                if Stepsize_options.speed_of_continuation
                    if n_rmv >= 3
                        Temp.speed_of_continuation = [];
                    else
                        Temp.speed_of_continuation = Temp.speed_of_continuation(1:end-n_rmv);
                    end
                end
                if Stepsize_options.rate_of_contraction
                    if n_rmv >= 3
                        Temp.rate_of_contraction = [];
                    else
                        Temp.rate_of_contraction = Temp.rate_of_contraction(1:end-n_rmv);
                    end
                end
                Plus = aux.clear_struct(Plus);
                Do.stepback = false;
                Do.suspend = false;
                Do.remove = true;
                Counter.remove = Counter.remove+1;
            else
                %% else
                Plus = aux.clear_struct(Plus);
                Do.stepback = false;
                Do.suspend = false;
                Do.remove = false;
            end
            if Opt.include_reverse && is_reverse && Solver.exitflag>0
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
                Jacobian.solver = aux.get_jacobian(fun,Path.var_all(:,end),Path.l_all(end),Opt);
            end
            [Bifurcation,Jacobian,Path] = bifurcation.check(fun,Jacobian,Path,Bifurcation,Info,res_corr,Solver,Opt);
        elseif aux.ison(Opt.bifurcation) && val && numel(Path.l_all)<=2
            if ~is_current_jacobian
                %% get jacobian if not current
                Jacobian.solver = aux.get_jacobian(fun,Path.var_all(:,end),Path.l_all(end),Opt);
            end
            Jacobian.sign_det = sign(det(Jacobian.solver(1:Info.nv,1:Info.nv)));
        end
        %
        %% step size control
        %
        % save latest stepsize
        %
        dsim1 = ds;
        %
        % save step size data
        %
        [Solver,Path] = aux.update_stepsize_data(Stepsize_options,Temp,Solver,Path);
        %
        % adjust stepsize
        %
        [ds,Counter,event_out] = step_size.control(ds,Counter,Solver,Do,Plus,Path,Jacobian,Opt,Info,event_out);
        %
        %% end loop
        %
        if val
            if ~isempty(Solver.output.iterations)
                aux.print_line(Opt,'-----> continued at l = %.4e\t|\tnew step size: ds = %.2e\t|\tloop counter = %d\t|\tstep = %d\t|\titerations = %d/%d\n',Path.l_all(end),ds,Counter.loop,Counter.step,Solver.output.iterations(end),Opt.n_iter_opt);
            else
                aux.print_line(Opt,'-----> continued at l = %.4e\t|\tnew step size: ds = %.2e\t|\tloop counter = %d\t|\tstep = %d\t|\titerations = %d/%d\n',Path.l_all(end),ds,Counter.loop,Counter.step,[],Opt.n_iter_opt);
            end
        else
            if ~isempty(Solver.output.iterations)
                aux.print_line(Opt,'-----> invalid point %s |\tnew step size: ds = %.2e\t|\tloop counter = %d\t|\tstep = %d\t|\titerations = %d/%d\n',inv_poi_str,ds,Counter.loop,Counter.step,Solver.output.iterations(end),Opt.n_iter_opt);
            else
                aux.print_line(Opt,'-----> invalid point %s |\tnew step size: ds = %.2e\t|\tloop counter = %d\t|\tstep = %d\t|\titerations = %d/%d\n',inv_poi_str,ds,Counter.loop,Counter.step,[],Opt.n_iter_opt);
            end
            
        end
        [Do,Info,Path,break_fun_out,Opt,Counter] = aux.exit_loop(Do,Info,Path,Opt,Counter,Bifurcation,ds,fun_solution,Jacobian,break_fun_out,val);
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
                [Plot, Opt] = plot.live_plot(Opt, Info, Path, ds, dsim1, Solver.output.iterations(end), Counter, fun_predictor, s_predictor, Plot, Bifurcation);
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
            Jacobian.last = [];
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
            plot.live_plot(Opt, Info, Path, ds, dsim1, Solver.output.iterations(end), Counter, fun_predictor, s_predictor, Plot, Bifurcation_last_plot);
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
    %% output
    %
    Info_out.number_of_steps = Counter.step;
    Info_out.number_of_invalid_points = Counter.loop - Counter.step;
    jacobian_out = Jacobian.last;
    exitflag = Info.exitflag;
    var_all = Path.var_all;
    l_all = Path.l_all;
    s_all = Path.s_all;
	%
end