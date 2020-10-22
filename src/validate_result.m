%% path continuation - validate_result
%
%   Institute of Dynamics and Vibration Research
%   Leibniz University Hannover
%   08.05.2020 - Alwin Förster
%
function [val,is_reverse] = validate_result(xs,fun_solution,vars,ls,ds,solver_exitflag,Opt)
    is_reverse = false;
    if solver_exitflag>0
        try
            if norm(fun_solution)<=Opt.solver_tol*10
                xi = [vars(:,end);ls(end)];
                norm_xs_xi = norm(xs-xi);
                if norm_xs_xi>=0.9*ds && norm_xs_xi<=1.1*ds
                    if length(ls)==1 && sign(xs(end)-ls)==sign(Opt.direction)
                        val = true;
                    else
                        xim1 = [vars(:,end-1);ls(end-1)];
                        alpha = acos(((xs-xi)'*(xi-xim1))/(sqrt((xs-xi)'*(xs-xi))*sqrt((xi-xim1)'*(xi-xim1))));
                        if alpha<Opt.alpha_reverse
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
        end
    else
        val = false;
    end
end