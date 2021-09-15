%% path continuation - validate_result
%
%   Institute of Dynamics and Vibration Research
%   Leibniz University Hannover
%   08.05.2020 - Alwin Förster
%
function [val,is_reverse,catch_flag] = validate_result(xs,fun_solution,var_all,l_all,ds,solver_exitflag,solver_jacobian,last_jacobian,do_convergeToTarget,bif,Opt)
    is_reverse = false;
    catch_flag = 0;
    if solver_exitflag>0
        try
            if ~Opt.check_residual || (norm(fun_solution)<=Opt.solver_tol*10)
                xi = [var_all(:,end);l_all(end)];
                norm_xs_xi = norm(xs-xi);
                if ((norm_xs_xi>=0.8*ds && norm_xs_xi<=1.2*ds || numel(l_all)==1)) || do_convergeToTarget || Opt.corrector.unique
                    if numel(l_all)==1
                        if numel(Opt.direction)==1 && sign(xs(end)-l_all(end))==sign(Opt.direction)
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
                        xim1 = [var_all(:,end-1);l_all(end-1)];
                        alpha = acos(((xs-xi)'*(xi-xim1))/(sqrt((xs-xi)'*(xs-xi))*sqrt((xi-xim1)'*(xi-xim1))));
                        if alpha<Opt.alpha_reverse
                            val = true;
                        else
                            if ~isempty(bif) && bif(1,end)==numel(l_all)
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