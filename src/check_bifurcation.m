%% path continuation - check_bifurcation
%
%   Institute of Dynamics and Vibration Research
%   Leibniz University Hannover
%   26.05.2020 - Alwin Förster
%
function [bif,sign_det_jacobian,bif_flag,var_all,l_all,s_all] = check_bifurcation(fun,solver_jacobian_red,var_all,l_all,s_all,bif,sign_det_jacobian,Opt)
    bif_flag = 0;
    solver_jacobian_red = full(solver_jacobian_red);
    if Opt.bifurcation.mark
        %% mark bifurcations:
        full_rank = length(solver_jacobian_red(:,1));
        rank_tol = 10^-2; % TODO!!!
        sign_det_current_jacobian = sign(det(solver_jacobian_red));
        if sign_det_current_jacobian*sign_det_jacobian<=0
            bif_type = (rank(solver_jacobian_red,rank_tol)==full_rank); % 1: true bif.; 0: zero point
            bif = [bif,[length(l_all);bif_type]];
            sign_det_jacobian = sign_det_current_jacobian;
            bif_flag = 1;
        end
    elseif Opt.bifurcation.determine || Opt.bifurcation.trace
        %% determine bifurcation-points:
        Opt_bif = Opt;
        Opt_bif.jacobian = false;
        [bif_solver,default_bif_solver_output] = continuation_solver(Opt_bif);
        det_solver_jacobian_red = det(solver_jacobian_red);
        residual_bif = @(x) residual_bifurcation(fun,x,Opt,1/det_solver_jacobian_red);
        full_rank = length(solver_jacobian_red(:,1));
        rank_tol = 10^-2; % TODO!!!
        sign_det_current_jacobian = sign(det_solver_jacobian_red);
        if sign_det_current_jacobian*sign_det_jacobian<=0
            %% find exact point:
            nds = 5;
            dss = linspace(s_all(end-1)-s_all(end),0,nds);
            dss = dss([3,2,4,1,5]);
            ind_bif = length(l_all);
            for i=1:nds
                dsp = dss(i);
                [var_bif_predictor,l_bif_predictor] = predictor(var_all,l_all,s_all,dsp,solver_jacobian_red,fun,Opt);
                dscale = get_dscale(Opt,var_bif_predictor,l_bif_predictor);
                [x_bif,fun_bif,bif_solver_exitflag,bif_solver_output,bif_solver_jacobian] = bif_solver(residual_bif,[var_bif_predictor;l_bif_predictor],dscale);
                if bif_solver_exitflag>0
                    s_all = [s_all(1:end-1),s_all(end-1)+[norm(x_bif-[var_all(:,end-1);l_all(end-1)]),norm(x_bif-[var_all(:,end-1);l_all(end-1)])+norm([var_all(:,end);l_all(end)]-x_bif)]];
                    var_all = [var_all(:,1:end-1),x_bif(1:end-1),var_all(:,end)];
                    l_all = [l_all(1:end-1),x_bif(end),l_all(end)];
                    solver_jacobian_red = bif_solver_jacobian;
                    full_rank = length(solver_jacobian_red(:,1));
                    break;
                end
            end
            bif_type = (rank(solver_jacobian_red,rank_tol)==full_rank); % 1: true bif.; 0: zero point
            bif = [bif,[ind_bif;bif_type]];
            sign_det_jacobian = sign_det_current_jacobian;
            bif_flag = 1;
        end
    else
        %% error
        error('No such bifurcation-mode');
    end
end