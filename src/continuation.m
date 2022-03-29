%% path continuation - continuation
%
%   Institute of Dynamics and Vibration Research
%   Leibniz University Hannover
%   08.05.2020 - Alwin F�rster
%
%   [var_all,l_all,exitflag,Bifurcation] = continuation(fun,var0,l_start,l_end,ds0,varargin)
%
%   fun = fun(var,l) != 0
%   l_start <= l <= l_end
%   ds0: initial stepsize
%
%% This file is part of continuation.
% 
% If you use continuation, please refer to:
%   A. F�rster, foerster@ids.uni-hannover.de
% 
% COPYRIGHT AND LICENSING: 
% Continuation Copyright (C) 2022 Alwin F�rster
%                                 (foerster@ids.uni-hannover.de)
%                                 Leibnitz University Hannover
% This program comes with NO WARRANTY. 
% Continuation is free software, you can redistribute and/or modify it
% under the GNU General Public License as published by the Free Software
% Foundation, either version 3 of the License, or (at your option) any
% later version. For details on license and warranty, see
% http://www.gnu.org/licenses or gpl-3.0.txt.
%
%%
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
    [Bifurcation,Counter,Do,Info,Info_out,Initial,Is,Jacobian,Path,Plot,Plus,Remove,Solver,Stepsize_information,Temp] = aux.initialize_structs(var0,l_start,l_end,ds0,Opt,Stepsize_options.rate_of_contraction);
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
        Is.current_jacobian = false;
        Counter.catch_old = Counter.catch;
        Temp = step_size.update_temp(Path,Solver,Stepsize_options,Temp);
        %
        %% residual
        %
        if Do.deflate
            try
                residual = @(x) aux.deflation(residual,x_deflation,x,Opt);
            catch
                aux.print_line(Opt,'---> delation: catch!\n');
                Counter.catch = Counter.catch + 1;
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
                x_predictor = homotopy.h_continuation(residual,[Path.var_all(:,end);Path.l_all(end)],Opt);
                var_predictor = x_predictor(1:end-1);
                l_predictor = x_predictor(end);
            else
                %% calc. predictor
                [var_predictor,l_predictor,fun_predictor,s_predictor,ds] = continuation.predictor(Path,ds,Jacobian.last,fun,res_corr,Solver,Opt);
                x_predictor = [var_predictor;l_predictor];
            end
        catch
            [var_predictor,l_predictor,fun_predictor,s_predictor,ds] = continuation.predictor(Path,ds,Jacobian.last,fun,res_corr,Solver,Opt);
            x_predictor = [var_predictor;l_predictor];
            aux.print_line(Opt,'---> predictor: catch!\n');
            Counter.catch = Counter.catch + 1;
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
            Is.current_jacobian = true;
        catch
            x_solution = NaN(size(x_predictor));
            fun_solution = inf(size(x_predictor));
            Solver.exitflag = -2;
            Solver.output = Solver.default_output;
            Do.convergeToTarget = false;
            aux.print_line(Opt,'---> solve: catch!\n');
            Counter.catch = Counter.catch + 1;
        end
        %
        %% calc stepsize information
        %
        % measure speed
        %
        if Stepsize_options.speed_of_continuation
            time_needed = toc(time_needed);
            Stepsize_information.speed_of_continuation = ds/time_needed; 
        end
        %
        % measure rate of contraction
        if Stepsize_options.rate_of_contraction
            if size(solver_stepsizes, 1) < 3
                Stepsize_information.rate_of_contraction = Opt.optimal_contraction_rate;
            else
                Stepsize_information.rate_of_contraction = solver_stepsizes(3,2)/solver_stepsizes(2,2);
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
        %
        %% adaptive corrector
        %
        [Do,Opt,corr_info] = corrector.adapt(Do,Opt,Path,Solver,fun,x_predictor,dscale,Jacobian.last,ds);
        %
        %% check result
        %
        % check result:
        [inv_poi_str,Do,Is,Opt] = ...
            aux.validate_result(x_solution,Plus,fun_solution,Path,ds,Solver,Jacobian,fun_predictor,s_predictor,Do,Bifurcation,Info,Is,Counter,Plot,Opt);
        % confirm result:
        [x_deflation,Bifurcation,Counter,Do,Info,Initial,Is,Jacobian,Path,Plus,Remove,Solver,Stepsize_information,Stepsize_options,Temp,Opt] = ...
            aux.confirm_result(fun,x_solution,x_predictor,Bifurcation,Counter,Do,Info,Initial,Is,Jacobian,Path,Plus,Remove,Solver,Stepsize_information,Stepsize_options,Temp,Opt);
        %
        %% Bifurcations
        %
        Bifurcation.flag = 0;
        if aux.ison(Opt.bifurcation) && Is.valid && ~Do.homotopy && numel(Path.l_all)>2
            if ~Is.current_jacobian
                %% get jacobian if not current
                Jacobian.solver = aux.get_jacobian(fun,Path.var_all(:,end),Path.l_all(end),Opt);
            end
            [Bifurcation,Jacobian,Path] = bifurcation.check(fun,Jacobian,Path,Bifurcation,Info,res_corr,Solver,Opt);
        elseif aux.ison(Opt.bifurcation) && Is.valid && numel(Path.l_all)<=2
            if ~Is.current_jacobian
                %% get jacobian if not current
                Jacobian.solver = aux.get_jacobian(fun,Path.var_all(:,end),Path.l_all(end),Opt);
            end
            Jacobian.sign_det = sign(det(Jacobian.solver(1:Info.nv,1:Info.nv)));
        end
        %
        %% step size control
        %
        % save latest stepsize:
        dsim1 = ds;
        % save step size data:
        [Solver,Path] = aux.update_stepsize_data(Stepsize_options,Temp,Solver,Path);
        % adjust stepsize:
        [ds,Counter,event_out,Opt] = step_size.control(ds,Counter,Solver,Do,Plus,Path,Jacobian,Opt,Info,event_out,Initial);
        %
        %% end loop
        %
        if Is.valid
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
        [Do,Info,Path,break_fun_out,Opt,Counter] = aux.exit_loop(Do,Info,Is,Path,Opt,Counter,Bifurcation,ds,fun_solution,Jacobian,break_fun_out);
        if Do.change_corrector
            Opt = aux.seton(Opt,'corrector',corr_info);
            res_corr = continuation.corrector(Opt);
            Do.change_corrector = false;
        end
        Opt = aux.update_Opt(Opt,Opt_is_set,Info);
        %
        %% live plot
        %
        if aux.ison(Opt.plot) && Is.valid
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