%% path continuation - corrector.adapt
%
%   Institute of Dynamics and Vibration Research
%   Leibniz University Hannover
%   20.01.2022 - Alwin Förster
%
function [Do,Opt,corr_info] = adapt(Do,Opt,Path,solver_exitflag,solver_output,Solver,fun,x_predictor,dscale,last_jacobian,ds)
    corr_info = [];
    if aux.ison(Opt.adapt_corrector)
        if Opt.corrector.sphere
            if solver_exitflag>0
                if Opt.adapt_corrector.solve
                    Opt_temp = aux.seton(Opt,'corrector','orthogonal');
                    res_corr_temp = continuation.corrector(Opt_temp);
                    residual_temp = @(x) aux.merge_residuals(fun,res_corr_temp,x,[Path.var_all;Path.l_all],ds,last_jacobian,Opt_temp);
                    [~,~,solver_exitflag_temp,solver_output_temp] = Solver.main(residual_temp,x_predictor,dscale);
                    if solver_exitflag_temp>0 && solver_output_temp.iterations<solver_output.iterations
                        corr_info = 'orthogonal';
                        Do.change_corrector = true;
                    end
                end
            else
                corr_info = 'orthogonal';
                Do.change_corrector = true;
            end
        elseif Opt.corrector.orthogonal
            if solver_exitflag>0
                if Opt.adapt_corrector.solve
                    Opt_temp = aux.seton(Opt,'corrector','sphere');
                    res_corr_temp = continuation.corrector(Opt_temp);
                    residual_temp = @(x) aux.merge_residuals(fun,res_corr_temp,x,[Path.var_all;Path.l_all],ds,last_jacobian,Opt_temp);
                    [~,~,solver_exitflag_temp,solver_output_temp] = Solver.main(residual_temp,x_predictor,dscale);
                    if solver_exitflag_temp>0 && solver_output_temp.iterations<solver_output.iterations
                        corr_info = 'sphere';
                        Do.change_corrector = true;
                    end
                end
            else
                corr_info = 'sphere';
                Do.change_corrector = true;
            end
        else
            corr_info = 'orthogonal';
            Do.change_corrector = true;
        end
    end
end