%% path continuation - predictor.adaptive
%
%   Institute of Dynamics and Vibration Research
%   Leibniz University Hannover
%   18.10.2020 - Alwin Förster
%
function [nt,nf] = adaptive(Path,Opt)
    if Opt.predictor_polynomial_adaptive
        var_old = Path.var_all(:,1:end-1);
        l_old = Path.l_all(1:end-1);
        s_old = Path.s_all(1:end-1);
        var_solution = Path.var_all(:,end);
        l_solution = Path.l_all(end);
        x_solution = [var_solution;l_solution];
        errmin = inf;
        ds_old = norm(x_solution-[var_old(:,end);l_old(end)]);
        for kt=1:Opt.predictor_polynomial_order
            for kf=0:Opt.predictor_polynomial_fit
                if length(l_old)==1
                    if numel(Opt.direction)==1
                        x_predictor_old = [var_old;l_old+sign(Opt.direction)*ds_old];
                    else
                        x_predictor_old = [var_old;l_old]+Opt.direction*ds_old;
                    end
                else
                    fpt = predictor.taylor(struct('var_all',var_old,'l_all',l_old,'s_all',s_old),kt,kf);
                    x_predictor_old = fpt(ds_old);
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
        nt = Opt.predictor_polynomial_order;
        nf = Opt.predictor_polynomial_fit;
    end
end