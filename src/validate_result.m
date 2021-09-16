%% path continuation - validate_result
%
%   Institute of Dynamics and Vibration Research
%   Leibniz University Hannover
%   08.05.2020 - Alwin Förster
%
function [val,is_reverse,catch_flag,Opt] = validate_result(x_solution,x_plus,fun_solution,Path,ds,solver_output,solver_exitflag,solver_jacobian,last_jacobian,fun_predictor,s_predictor,Do,Bifurcation,Info,Counter,Plot,Opt)
    %% automated validation
    %
    is_reverse = false;
    catch_flag = 0;
    if solver_exitflag>0
        try
            if ~Opt.check_residual || (norm(fun_solution)<=Opt.solver_tol*10)
                xi = [Path.var_all(:,end);Path.l_all(end)];
                norm_xs_xi = norm(x_solution-xi);
                if ((norm_xs_xi>=0.8*ds && norm_xs_xi<=1.2*ds || numel(Path.l_all)==1)) || Do.convergeToTarget || Opt.corrector.unique
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
                            end
                        end
                    end
                else
                    val = false;
                end
            else
                val = false;
            end
        catch
            val = false;
            catch_flag = 1;
        end
    else
        val = false;
    end
    %
    %% approve manually
    if val && ((islogical(Opt.approve_manually) && Opt.approve_manually) || (~islogical(Opt.approve_manually) && (x_solution(end)-Opt.approve_manually)*(Path.l_all(end)-Opt.approve_manually)<=0))
        Opt.approve_manually = true;
        if numel(Path.l_all)>1
            try
                Path_app = Path;
                if isempty(x_plus)
                    Path_app.var_all = [Path.var_all,x_solution(1:end-1)];
                    Path_app.l_all = [Path.l_all,x_solution(end)];
                    Path_app.s_all = [Path.s_all,Path.s_all(end)+norm(x_solution-[Path.var_all(:,end-1);Path.l_all(end-1)])];
                else
                    Path_app.var_all = [Path.var_all,x_solution(1:end-1),x_plus(1:end-1)];
                    Path_app.l_all = [Path.l_all,x_solution(end),x_plus(end)];
                    Path_app.s_all = [Path.s_all,Path.s_all(end)+norm(x_solution-[Path.var_all(:,end-2);Path.l_all(end-2)])*[1,1]+norm(x_plus-x_solution)*[0,1]];
                end
                [Plot, Opt] = live_plot(Opt, Info, Path, ds, ds, solver_output.iterations, Counter, fun_predictor, s_predictor, Plot, Bifurcation);
            catch
                warning('The plot update for approval has failed.');
            end
        end
        prompt = sprintf('------> approve point at l = %.4e (y/n): ',Path.l_all(end));
        correct_input = false;
        while ~correct_input
            y_or_n = input(prompt,'s');
            if strcmp(y_or_n,'y')
                correct_input = true;
            elseif strcmp(y_or_n,'n')
                correct_input = true;
                val = false;
            elseif strcmp(y_or_n,'off')
                correct_input = true;
                Opt.approve_manually = false;
            elseif ~isnan(str2double(y_or_n))
                correct_input = true;
                Opt.approve_manually = str2double(y_or_n);
            else
                fprintf('------> Enter ''y'' for yes or ''n'' for no! (Deactivate with ''off'' or set double limit)\n');
            end
        end
    end
    %
end