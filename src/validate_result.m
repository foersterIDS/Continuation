%% path continuation - validate_result
%
%   Institute of Dynamics and Vibration Research
%   Leibniz University Hannover
%   08.05.2020 - Alwin Förster
%
function [val,is_reverse,catch_flag] = validate_result(xs,fun_solution,vars,ls,ds,solver_exitflag,solver_jacobian,last_jacobian,do_convergeToTarget,Opt)
    is_reverse = false;
    catch_flag = 0;
    if solver_exitflag>0
        try
            if norm(fun_solution)<=Opt.solver_tol*10
                xi = [vars(:,end);ls(end)];
                norm_xs_xi = norm(xs-xi);
                if ((norm_xs_xi>=0.8*ds && norm_xs_xi<=1.2*ds || numel(ls)==1)) || do_convergeToTarget || Opt.corrector.unique
                    if length(ls)==1
                        if numel(Opt.direction)==1 && sign(xs(end)-ls)==sign(Opt.direction)
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
                        ori_reverse = true;
                        x_all = [vars; ls];
                        if numel(ls) > 3 && sum(abs(size(solver_jacobian) - size(last_jacobian)))==0
                            if diff(size(solver_jacobian)) == 0
                                ori_new = sign(det([solver_jacobian(1:end-1,:).', x_all(:,end)-x_all(:,end-1)])); % solver_jacobian
                            else
                                ori_new = sign(det([solver_jacobian.', x_all(:,end)-x_all(:,end-1)])); % solver_jacobian
                            end
                            if diff(size(last_jacobian)) == 0
                                ori_old = sign(det([last_jacobian(1:end-1,:).', x_all(:,end-1)-x_all(:,end-2)])); % last_jacobian
                            else
                                ori_old = sign(det([last_jacobian.', x_all(:,end-1)-x_all(:,end-2)])); % last_jacobian
                            end
                            if ori_old*ori_new < 0
                                ori_reverse = false;
                            end
                        end
                        xim1 = [vars(:,end-1);ls(end-1)];
                        alpha = acos(((xs-xi)'*(xi-xim1))/(sqrt((xs-xi)'*(xs-xi))*sqrt((xi-xim1)'*(xi-xim1))));
                        if alpha<Opt.alpha_reverse && ori_reverse
                            val = true;
                        else
                            val = false;
                            is_reverse = true;
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