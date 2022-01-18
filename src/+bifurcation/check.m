%% path continuation - bifurcation.check
%
%   Institute of Dynamics and Vibration Research
%   Leibniz University Hannover
%   26.05.2020 - Alwin Förster
%   02.07.2021 - Tido Kubatschek
%
function [Bifurcation,sign_det_jacobian,Path] = check(fun,solver_jacobian_red,Path,Bifurcation,sign_det_jacobian,res_corr,Solver,Opt)
    Bifurcation.flag = 0;
    solver_jacobian_red = full(solver_jacobian_red);
    if Opt.bifurcation.mark
        %% mark bifurcations:
        sign_det_current_jacobian = sign(det(solver_jacobian_red));
        if sign_det_current_jacobian*sign_det_jacobian<=0
            bif_type = NaN; % 1: fold bif.; 0: branch point bif; NaN: unknown
            Bifurcation.bif = [Bifurcation.bif,[numel(Path.l_all);bif_type]];
            sign_det_jacobian = sign_det_current_jacobian;
            Bifurcation.flag = 1;
        end
    elseif Opt.bifurcation.determine || Opt.bifurcation.trace
        %% determine bifurcation-points:
        Opt_bif = Opt;
        Opt_bif.jacobian = false;
        [bif_solver,default_bif_solver_output] = continuation.solver(Opt_bif);
        det_solver_jacobian_red = det(solver_jacobian_red);
        residual_bif = @(x) bifurcation.residual(fun,x,Opt,1/det_solver_jacobian_red);
        full_rank = length(solver_jacobian_red(:,1));
        rank_tol = Opt_bif.solver_tol * 10000; % TODO!!!
        sign_det_current_jacobian = sign(det_solver_jacobian_red);
        if sign_det_current_jacobian*sign_det_jacobian<=0
            %% find exact point:
            nds = 5;
            dss = linspace(Path.s_all(end-1)-Path.s_all(end),0,nds);
            dss = dss([3,2,4,1,5]);
            ind_bif = length(Path.l_all);
            bif_type = NaN;
            for i=1:nds
                dsp = dss(i);
                [var_bif_predictor,l_bif_predictor] = continuation.predictor(Path,dsp,solver_jacobian_red,fun,res_corr,Solver,Opt);
                dscale = aux.get_dscale(Opt,struct('var_all',var_bif_predictor,'l_all',l_bif_predictor));
                [x_bif,fun_bif,bif_solver_exitflag,bif_solver_output,bif_solver_jacobian] = bif_solver(residual_bif,[var_bif_predictor;l_bif_predictor],dscale);
                if bif_solver_exitflag>0
                    Path.s_all = [Path.s_all(1:end-1),Path.s_all(end-1)+[norm(x_bif-[Path.var_all(:,end-1);Path.l_all(end-1)]),norm(x_bif-[Path.var_all(:,end-1);Path.l_all(end-1)])+norm([Path.var_all(:,end);Path.l_all(end)]-x_bif)]];
                    Path.var_all = [Path.var_all(:,1:end-1),x_bif(1:end-1),Path.var_all(:,end)];
                    Path.l_all = [Path.l_all(1:end-1),x_bif(end),Path.l_all(end)];
                    nv = numel(Path.var_all(:,1));
                    solver_jacobian_red = bif_solver_jacobian(1:nv,1:nv);
                    solver_jacobian_lam = bif_solver_jacobian(1:nv,nv+1);
                    full_rank = length(solver_jacobian_red(:,1));
                    %
                    % get type of bifurcation
                    jac_red_jac_lam = [solver_jacobian_red, solver_jacobian_lam];
                    bif_type = (rank(jac_red_jac_lam,rank_tol)==full_rank); % 1: fold bif.; 0: branch point bif; NaN: unknown
                    %
                    % if bif_type = 0 (branch point) calculate directions
                    % of paths by null() and save to Bifurcation.dirs cell array
                    if bif_type == 0
                        null_matrix = null(jac_red_jac_lam);
                        for k_dirs = 1:width(null_matrix)
                            bif_dir = null_matrix(:,k_dirs);
                            if height(Bifurcation.dirs) == 1 && isempty(Bifurcation.dirs{1,2})
                                Bifurcation.dirs(1,:) = [{ind_bif}, {bif_dir}];
                            else
                                Bifurcation.dirs(end+1,:) = [{ind_bif}, {bif_dir}];
                            end
                        end
                    end
                    break;
                end
            end
            Bifurcation.bif = [Bifurcation.bif,[ind_bif;bif_type]];
            sign_det_jacobian = sign_det_current_jacobian;
            Bifurcation.flag = 1;
        end
    else
        %% error
        error('No such bifurcation-mode');
    end
end