%% path continuation - validate_result
%
%   Institute of Dynamics and Vibration Research
%   Leibniz University Hannover
%   08.05.2020 - Alwin Förster
%
function [val,is_reverse] = validate_result(xs,vars,ls,solver_exitflag,Opt)
    is_reverse = false;
    if solver_exitflag>0
        try
            if length(ls)==1 && sign(xs(end)-ls)==sign(Opt.direction)
                val = true;
            else
                xi = [vars(:,end);ls(end)];
                xim1 = [vars(:,end-1);ls(end-1)];
                alpha = acos(((xs-xi)'*(xi-xim1))/(sqrt((xs-xi)'*(xs-xi))*sqrt((xi-xim1)'*(xi-xim1))));
                if alpha<Opt.alpha_reverse
                    val = true;
                else
                    val = false;
                    is_reverse = true;
                end
            end
        catch
            val = false;
        end
    else
        val = false;
    end
end