%% path continuation - aux.validate_result
%
%   Institute of Dynamics and Vibration Research
%   Leibniz University Hannover
%   08.05.2020 - Alwin Förster
%
function [val,is_reverse,catch_flag,inv_poi_str,Do,Opt] = validate_result(x_solution,Plus,fun_solution,Path,ds,solver_output,solver_exitflag,Jacobian,fun_predictor,s_predictor,Do,Bifurcation,Info,Counter,Plot,Opt)
    %% automated validation
    %
    is_reverse = false;
    catch_flag = 0;
    inv_poi_str = '                          ';
    if solver_exitflag>0
        try
            if ~Opt.check_residual || (norm(fun_solution)<=Opt.solver_tol*10)
                xi = [Path.var_all(:,end);Path.l_all(end)];
                norm_xs_xi = norm(x_solution-xi);
                if ((norm_xs_xi>=Opt.ds_tol(1)*ds && norm_xs_xi<=Opt.ds_tol(2)*ds || numel(Path.l_all)==1) && norm_xs_xi<=Opt.ds_tol(2)*Opt.ds_max) || Do.convergeToTarget || Opt.corrector.unique
                    if numel(Path.l_all)==1
                        if numel(Opt.direction)==1 && sign(x_solution(end)-Path.l_all(end))==sign(Opt.direction)
                            val = true;
                        else
                            alpha = acos(((x_solution-xi)'*Opt.direction)/(sqrt((x_solution-xi)'*(x_solution-xi))*sqrt(Opt.direction'*Opt.direction)));
                            if alpha<Opt.alpha_reverse
                                val = true;
                            else
                                val = false;
                                is_reverse = true;
                                inv_poi_str_temp = sprintf('(alpha: %.2f)',alpha);
                                inv_poi_str(1:numel(inv_poi_str_temp)) = inv_poi_str_temp;
                            end
                        end
                    else
                        xim1 = [Path.var_all(:,end-1);Path.l_all(end-1)];
                        alpha = acos(((x_solution-xi)'*(xi-xim1))/(sqrt((x_solution-xi)'*(x_solution-xi))*sqrt((xi-xim1)'*(xi-xim1))));
                        if alpha<Opt.alpha_reverse
                            val = true;
                        else
                            if ~isempty(Bifurcation.bif) && Bifurcation.bif(1,end)==numel(Path.l_all)
                                val = true;
                            else
                                val = false;
                                is_reverse = true;
                                inv_poi_str_temp = sprintf('(alpha: %.2f)',alpha);
                                inv_poi_str(1:numel(inv_poi_str_temp)) = inv_poi_str_temp;
                            end
                        end
                    end
                else
                    val = false;
                    inv_poi_str(1:25) = '(dx outside of ds-bounds)';
                end
            else
                val = false;
                inv_poi_str(1:18) = '(fun=0 not solved)';
            end
        catch
            val = false;
            catch_flag = 1;
            inv_poi_str(1:24) = '(error while validating)';
        end
    else
        val = false;
        inv_poi_str(1:26) = '(Solver does not converge)';
    end
    %
    %% enforce_ds_max
    %
    if val && Opt.enforce_ds_max
        val = logical(heaviside(min(Opt.ds_max-(x_solution-xi))));
        if ~val
            inv_poi_str(1:17) = '(ds_max exceeded)';
        end
    end
    %
    %% approve manually
    %
    if val && ((islogical(Opt.approve_manually) && Opt.approve_manually) || (~islogical(Opt.approve_manually) && (x_solution(end)-Opt.approve_manually)*(Path.l_all(end)-Opt.approve_manually)<=0))
        Opt.approve_manually = true;
        if numel(Path.l_all)>1
            try
                Path_app = Path;
                if isempty(Plus.x)
                    Path_app.var_all = [Path.var_all,x_solution(1:end-1)];
                    Path_app.l_all = [Path.l_all,x_solution(end)];
                    Path_app.s_all = [Path.s_all,Path.s_all(end)+norm(x_solution-[Path.var_all(:,end-1);Path.l_all(end-1)])];
                else
                    Path_app.var_all = [Path.var_all,x_solution(1:end-1),Plus.x(1:end-1)];
                    Path_app.l_all = [Path.l_all,x_solution(end),Plus.x(end)];
                    Path_app.s_all = [Path.s_all,Path.s_all(end)+norm(x_solution-[Path.var_all(:,end-2);Path.l_all(end-2)])*[1,1]+norm(Plus.x-x_solution)*[0,1]];
                end
                [Plot, Opt] = live_plot(Opt, Info, Path, ds, ds, solver_output.iterations, Counter, fun_predictor, s_predictor, Plot, Bifurcation);
            catch
                aux.print_line(Opt,'--> The plot update for approval has failed.\n');
            end
        end
        prompt = sprintf('------> approve point at l = %.4e (y/n): ',Path.l_all(end));
        correct_input = false;
        while ~correct_input
            input_string = input(prompt,'s');
            if strcmp(input_string,'y')
                correct_input = true;
            elseif strcmp(input_string,'n')
                correct_input = true;
                val = false;
                inv_poi_str(1:19) = '(rejected manually)';
            elseif strcmp(input_string,'off')
                correct_input = true;
                Opt.approve_manually = false;
            elseif strcmp(input_string,'exit')
                correct_input = true;
                val = false;
                Do.stop_manually = true;
                inv_poi_str(1:19) = '(rejected manually)';
            elseif ~isnan(str2double(input_string))
                correct_input = true;
                Opt.approve_manually = str2double(input_string);
            else
                aux.print_line(Opt,'------> Enter ''y'' for yes or ''n'' for no! (Deactivate with ''off'', set double limit or leave using ''exit'')\n');
            end
        end
    end
    %
end