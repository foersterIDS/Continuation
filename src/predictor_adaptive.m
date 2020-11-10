%% path continuation - predictor_adaptive
%
%   Institute of Dynamics and Vibration Research
%   Leibniz University Hannover
%   18.10.2020 - Alwin Förster
%
function [nt,nf] = predictor_adaptive(var_all,l_all,s_all,Opt)
    if Opt.predictor_adaptive
        var_old = var_all(:,1:end-1);
        l_old = l_all(1:end-1);
        s_old = s_all(1:end-1);
        var_solution = var_all(:,end);
        l_solution = l_all(end);
        x_solution = [var_solution;l_solution];
        errmin = inf;
        ds_old = norm(x_solution-[var_old(:,end);l_old(end)]);
        for kt=1:Opt.predictor_taylor
            for kf=0:Opt.predictor_fit
                if length(l_old)==1
                    if numel(Opt.direction)==1
                        x_predictor_old = [var_old;l_old+sign(Opt.direction)*ds_old];
                    else
                        x_predictor_old = [var_old;l_old]+Opt.direction*ds_old;
                    end
                else
                    x_predictor_old = predictor_taylor(var_old,l_old,s_old,kt,kf,ds_old);
                end
                err = norm(x_predictor_old-x_solution);
                if err<errmin
                    errmin = err;
                    nt = kt;
                    nf = kf;
                end
            end
        end
    else
        nt = Opt.predictor_taylor;
        nf = Opt.predictor_fit;
    end
end