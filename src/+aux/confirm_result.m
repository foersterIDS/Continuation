%% path continuation - aux.confirm_result
%
%   Institute of Dynamics and Vibration Research
%   Leibniz University Hannover
%   28.03.2022 - Alwin Förster
%   12.08.2022 - Anna Lefken
%
function [x_deflation,Bifurcation,Counter,Do,Info,Initial,Is,Jacobian,Path,Plus,Remove,Solver,Stepsize_information,Stepsize_options,Temp,Opt] = ...
    confirm_result(func,x_solution,x_predictor,Bifurcation,Counter,Do,Info,Initial,Is,Jacobian,Path,Plus,Remove,Solver,Stepsize_information,Stepsize_options,Temp,Opt,Opt_is_set)
    x_deflation = x_solution;
    if Is.valid
        %% valid result
        if isempty(Plus.x)
            Path.var_all = [Path.var_all,x_solution(1:end-1)];
            Path.l_all = [Path.l_all,x_solution(end)];
            Path.s_all = [Path.s_all,Path.s_all(end)+norm(x_solution-[Path.var_all(:,end-1);Path.l_all(end-1)])];
            if Opt_is_set.bif_additional_testfunction
                Path.biftest_value=[Path.biftest_value Opt.bif_additional_testfunction(func,x_solution,Jacobian,Path,Info)];
            end
            % update stepsize information
            % measure speed
            %
            if Stepsize_options.speed_of_continuation
                if ~isempty(Temp.speed_of_continuation)
                    Temp.speed_of_continuation = [Temp.speed_of_continuation, Stepsize_information.speed_of_continuation];
                else
                    Temp.speed_of_continuation = Stepsize_information.speed_of_continuation;
                end
            end
            %
            % measure rate of contraction
            if Stepsize_options.rate_of_contraction
                if ~isempty(Temp.rate_of_contraction)
                    Temp.rate_of_contraction = [Temp.rate_of_contraction, Stepsize_information.rate_of_contraction];
                else
                    Temp.rate_of_contraction = Stepsize_information.rate_of_contraction;
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
            if Opt_is_set.bif_additional_testfunction
                Path.biftest_value=[Path.biftest_value,Opt.bif_additional_testfunction(func,x_solution,Jacobian,Path,Info),Plus.biftest_value];
            end
            
            if Stepsize_options.predictor
                Temp.predictor = [Temp.predictor, x_predictor, Plus.x_predictor];
                Plus.x_predictor = [];
            end
            if Stepsize_options.speed_of_continuation
                Temp.speed_of_continuation = [Temp.speed_of_continuation, Stepsize_information.speed_of_continuation, Plus.speed_of_continuation];
                Plus.speed_of_continuation = [];
            end
            if Stepsize_options.rate_of_contraction
                Temp.rate_of_contraction = [Temp.rate_of_contraction, Stepsize_information.rate_of_contraction, Plus.rate_of_contraction];
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
            Jacobian.last = aux.get_jacobian(func,Path.var_all(:,end),Path.l_all(end),Opt);
        end
        if Do.remove
            if Path.s_all(end)>(Remove.s+Remove.ds)
                Opt.ds_max = Initial.ds_max;
                Do.remove = false;
                Counter.remove = 0;
            end
        end
    else
        %% invalid result
        Counter.error = Counter.error+1;
        if Opt.deflation && ~Do.deflate && ~isnan(sum(x_solution(:,end))) && Is.reverse && Counter.error>=Opt.deflation_error_counter
            Do.deflate = true;
            x_deflation = x_solution;
        else
            Do.deflate = false;
        end
        if ((Counter.error==Opt.stepback_error_counter) || Do.stepback_manually) && (numel(Path.l_all)>1)
            %% stepback
            if Counter.valid_stepback<Opt.stepback_error_counter
                Do.stepback = true;
                Plus.x = [Path.var_all(:,end);Path.l_all(end)];
                if Opt_is_set.bif_additional_testfunction
                    Plus.biftest_value=Opt.bif_additional_testfunction(func,x_solution,Jacobian,Path,Info);
                end
                if Stepsize_options.predictor
                    Plus.x_predictor = x_predictor;
                end
                if Stepsize_options.speed_of_continuation
                    Plus.speed_of_continuation = Stepsize_information.speed_of_continuation;
                end
                if Stepsize_options.rate_of_contraction
                    Plus.rate_of_contraction = Stepsize_information.rate_of_contraction;
                end
            else
                Do.stepback = false;
            end
            Path.var_all(:,end) = [];
            Path.l_all(end) = [];
            Path.s_all(end) = [];
            if Opt_is_set.bif_additional_testfunction
                Path.biftest_value(end) = [];
            end
            if Stepsize_options.predictor
                if ~isempty(Temp.predictor)
                    Temp.predictor(:,end) = [];
                else
                    Temp.predictor = [];
                end
            end
            if Stepsize_options.speed_of_continuation && ~isempty(Temp.speed_of_continuation)
                Temp.speed_of_continuation(end) = [];
            end
            if Stepsize_options.rate_of_contraction && ~isempty(Temp.rate_of_contraction)
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
                if Opt_is_set.bif_additional_testfunction
                    Path.biftest_value = [Path.biftest_value,Plus.biftest_value];
                end
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
            Remove.s = Path.s_all(n_path);
            n_rmv = min([2*Opt.remove_error_counter,n_path-1]);
            Opt.ds_max = max([mean(diff(Path.s_all(n_path+((-ceil(n_rmv/2)+1):0))))*0.75,Opt.ds_min]);
            Path.var_all(:,n_path+((-n_rmv+1):0)) = [];
            Path.l_all(n_path+((-n_rmv+1):0)) = [];
            Path.s_all(n_path+((-n_rmv+1):0)) = [];
            if Opt_is_set.bif_additional_testfunction
                Path.biftest_value(n_path+((-n_rmv+1):0)) = [];
            end
            Remove.ds = Remove.s-Path.s_all(end);
            dscale_rmv = aux.get_dscale(Opt,Path);
            residual_fixed_value_rmv = @(v) aux.residual_fixed_value(func,v,Path.l_all(end),Opt);
            [var_rmv,fun_rmv,~,~,rmv_jacobian] = Solver.main(residual_fixed_value_rmv,Path.var_all(:,end),dscale_rmv(1:end-1));
            if ~isempty(Bifurcation.bif) && numel(Bifurcation.bif(1,:))>0
                n_bifs_rmv = sum(sum(Bifurcation.bif(1,:)'==(n_path+((-n_rmv+1):0))));
                if n_bifs_rmv>0
                    Bifurcation.bif(:,end+((1-n_bifs_rmv):0)) = [];
                    Jacobian.sign_det_red = Jacobian.sign_det_red*(-1)^(n_bifs_rmv);
                end
            end
            Path.var_all(:,end) = var_rmv;
            Jacobian.last = [rmv_jacobian,aux.numeric_jacobian(@(x) func(x(1:Info.nv),x(Info.nv+1)),[Path.var_all(:,end);Path.l_all(end)],'central_value',fun_rmv,'derivative_dimensions',Info.nv+1,'diffquot',Opt.diffquot)];
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
        if Opt.include_reverse && Is.reverse && Solver.exitflag>0
            Path = aux.include_reverse(x_solution,Path);
        end
        if aux.ison(Opt.homotopy) && Counter.error>=Opt.homotopy_error_counter
            Do.homotopy = true;
        else
            Do.homotopy = false;
        end
        if Is.catch
            Counter.catch = Counter.catch + 1;
        end
    end
end