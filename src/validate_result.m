%% path continuation - validate_result
%
%   Institute of Dynamics and Vibration Research
%   Leibniz University Hannover
%   08.05.2020 - Alwin Förster
%
function [val,is_reverse,catch_flag] = validate_result(xs,fun_solution,Path,ds,solver_exitflag,solver_jacobian,last_jacobian,Do,Bifurcation,Opt)
    is_reverse = false;
    catch_flag = 0;
    if solver_exitflag>0
        try
            if ~Opt.check_residual || (norm(fun_solution)<=Opt.solver_tol*10)
                xi = [Path.var_all(:,end);Path.l_all(end)];
                norm_xs_xi = norm(xs-xi);
                if ((norm_xs_xi>=0.8*ds && norm_xs_xi<=1.2*ds || numel(Path.l_all)==1)) || Do.convergeToTarget || Opt.corrector.unique
                    if numel(Path.l_all)==1
                        if numel(Opt.direction)==1 && sign(xs(end)-Path.l_all(end))==sign(Opt.direction)
                            val = true;
                        else
                            alpha = acos(((xs-xi)'*Opt.direction)/(sqrt((xs-xi)'*(xs-xi))*sqrt(Opt.direction'*Opt.direction)));
                            if alpha<Opt.alpha_reverse
                                val = true;
                            else
                                val = false;
                                is_reverse = true;
                            end
                        end
                    else
                        xim1 = [Path.var_all(:,end-1);Path.l_all(end-1)];
                        alpha = acos(((xs-xi)'*(xi-xim1))/(sqrt((xs-xi)'*(xs-xi))*sqrt((xi-xim1)'*(xi-xim1))));
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
end