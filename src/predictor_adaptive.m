%% path continuation - predictor_adaptive
%
%   Institute of Dynamics and Vibration Research
%   Leibniz University Hannover
%   18.10.2020 - Alwin Förster
%
function [nt,nf] = predictor_adaptive(vars,ls,Opt)
    if Opt.predictor_adaptive
        vars_old = vars(:,1:end-1);
        ls_old = ls(1:end-1);
        var_solution = vars(:,end);
        l_solution = ls(end);
        x_solution = [var_solution;l_solution];
        errmin = inf;
        ds_old = norm(x_solution-[vars_old(:,end);ls_old(end)]);
        for kt=1:Opt.predictor_taylor
            for kf=0:Opt.predictor_fit
                x_predictor_old = predictor_taylor(vars_old,ls_old,kt,kf,ds_old);
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