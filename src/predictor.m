%% path continuation - predictor
%
%   Institute of Dynamics and Vibration Research
%   Leibniz University Hannover
%   08.05.2020 - Alwin Förster
%
function [vp,lp] = predictor(var_all,l_all,s_all,ds,solver_jacobian,fun,Opt)
    if Opt.predictor.polynomial
        if length(l_all)==1
            if numel(Opt.direction)==1
                xip1 = [var_all;l_all+sign(Opt.direction)*ds];
            else
                xip1 = [var_all;l_all]+Opt.direction*ds;
            end
        else
            [nt,nf] = predictor_adaptive(var_all,l_all,s_all,Opt);
            xip1 = predictor_taylor(var_all,l_all,s_all,nt,nf,ds);
        end
    elseif Opt.predictor.tangential
        if length(l_all)==1
            if numel(Opt.direction)==1
                xip1 = [var_all;l_all+sign(Opt.direction)*ds];
            else
                xip1 = [var_all;l_all]+Opt.direction*ds;
            end
        else
            xip1 = predictor_ode(var_all,l_all,ds,solver_jacobian,fun);
        end
    else
        error('predictor not set or of unknown type');
    end
    vp = xip1(1:end-1);
    lp = xip1(end);
end