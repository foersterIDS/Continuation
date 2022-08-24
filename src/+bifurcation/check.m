%% path continuation - bifurcation.check
%
%   Institute of Dynamics and Vibration Research
%   Leibniz University Hannover
%   26.05.2020 - Alwin Förster
%   02.07.2021 - Tido Kubatschek
%
function [Bifurcation,Jacobian,Path] = check(func,Jacobian,Path,Bifurcation,Info,res_corr,Solver,Opt,Opt_is_set)
    Bifurcation.flag = 0;
    last_jacobian_red = Jacobian.last(1:Info.nv,1:Info.nv);
    last_jacobian_red = full(last_jacobian_red);
    bif_found=0;    % skip additional testfunction, when det(jac)=0
    if Opt.bifurcation.mark
        %% mark bifurcations:
        sign_det_current_jacobian_red = sign(det(last_jacobian_red));
        if sign_det_current_jacobian_red*Jacobian.sign_det_red<=0
            bif_found=1;
            sign_det_current_jacobian = sign(det(Jacobian.last));
            bif_type = (sign(det(Jacobian.previous))==sign(det(Jacobian.last))); % 1: fold bif.; 0: branch point bif; NaN: unknown
            Bifurcation.bif = [Bifurcation.bif,[numel(Path.l_all);bif_type]];
            Jacobian.sign_det_red = sign_det_current_jacobian_red;
            Jacobian.sign_det = sign_det_current_jacobian;
            Bifurcation.flag = 1;
            Bifurcation.scaling = [Bifurcation.scaling,1];
        end
    elseif Opt.bifurcation.determine || Opt.bifurcation.trace || Opt.bifurcation.parameter_trace
        %% determine bifurcation-points:
        Opt_bif = Opt;
        Opt_bif.jacobian = false;
        [bif_solver,default_bif_solver_output] = continuation.solver(Opt_bif,0);
        det_solver_jacobian_red = det(last_jacobian_red);
        Bifurcation.scaling = [Bifurcation.scaling,1/det_solver_jacobian_red];
        residual_bif = @(x) bifurcation.residual(func,x,Opt,Info,Bifurcation.scaling(end));
        sign_det_current_jacobian_red = sign(det_solver_jacobian_red);
        if sign_det_current_jacobian_red*Jacobian.sign_det_red<=0
            bif_found=1;
            %% find exact point:
            nds = 1000;
            dss = linspace(Path.s_all(end-1)-Path.s_all(end),0,nds);
            ind_bif = length(Path.l_all);
            bif_type = NaN;
            nv = numel(Path.var_all(:,1));
            
            det_left=det(Jacobian.previous(1:nv,1:nv));
            det_right=det(Jacobian.last(1:nv,1:nv));
            start_dsp=det_right/(det_right-det_left)*(dss(1)-dss(end));     % estimated bifurcation point
            ind_start=find(dss>start_dsp,1);                                % index of estimated bifurcation point
            
            ind_dss=zeros(nds,1);                                           % index vector for sorting dss
            ind_dss(1)=ind_start;                                           % alternating left and right of ind_start
            if ind_start<nds/2
                ind_alternating_less=ind_start-1:-1:1;
                ind_dss(2:2:length(ind_alternating_less)*2)=ind_alternating_less;
                ind_dss(3:2:length(ind_alternating_less)*2+1)=ind_start+1:1:ind_start+length(ind_alternating_less);
                ind_dss(length(ind_alternating_less)*2+2:end)=ind_start+length(ind_alternating_less)+1:1:nds;
            else
                ind_alternating_greater=ind_start+1:1:nds;
                ind_dss(2:2:length(ind_alternating_greater)*2)=ind_start-1:-1:ind_start-length(ind_alternating_greater);
                ind_dss(3:2:length(ind_alternating_greater)*2+1)=ind_alternating_greater;
                ind_dss(length(ind_alternating_greater)*2+2:nds)=ind_start-length(ind_alternating_greater)-1:-1:1;
            end
            dss=dss(ind_dss);                                               % sort dss
            
            for i=1:nds
                dsp = dss(i);
                [var_bif_predictor,l_bif_predictor] = continuation.predictor(Path,dsp,last_jacobian_red,func,res_corr,Solver,Opt);
                dscale = aux.get_dscale(Opt,struct('var_all',var_bif_predictor,'l_all',l_bif_predictor));
                [x_bif,fun_bif,bif_solver_exitflag,bif_solver_output,bif_solver_jacobian] = bif_solver(residual_bif,[var_bif_predictor;l_bif_predictor],dscale);
                if bif_solver_exitflag>0
                    Path.s_all = [Path.s_all(1:end-1),Path.s_all(end-1)+[norm(x_bif-[Path.var_all(:,end-1);Path.l_all(end-1)]),norm(x_bif-[Path.var_all(:,end-1);Path.l_all(end-1)])+norm([Path.var_all(:,end);Path.l_all(end)]-x_bif)]];
                    Path.var_all = [Path.var_all(:,1:end-1),x_bif(1:end-1),Path.var_all(:,end)];
                    Path.l_all = [Path.l_all(1:end-1),x_bif(end),Path.l_all(end)];
                    if Opt_is_set.bif_additional_testfunction
                        Path.biftest_value=[Path.biftest_value Opt.bif_additional_testfunction(func,x_bif,Jacobian,Path,Info)];
                    end
                    % 
                    % get type of bifurcation
                    bif_type = (sign(det(Jacobian.previous))==sign(det(Jacobian.last))); % 1: fold bif.; 0: branch point bif; NaN: unknown
                    %
                    % if bif_type = 0 (branch point) calculate directions
                    % of paths by null() and save to Bifurcation.dirs cell array
                    if bif_type == 0
                        nv = numel(Path.var_all(:,1));
                        solver_jacobian_red = bif_solver_jacobian(1:nv,1:nv);
                        solver_jacobian_lam = bif_solver_jacobian(1:nv,nv+1);
                        jac_red_jac_lam = [solver_jacobian_red, solver_jacobian_lam];
                        %
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
            Jacobian.sign_det_red = sign_det_current_jacobian_red;
            Bifurcation.flag = 1;
        end
    else
        %% error
        error('No such bifurcation-mode');
    end
    if Opt_is_set.bif_additional_testfunction && bif_found==0
        if Opt.bifurcation.mark
            %% mark bifurcations:
            if Path.biftest_value(end-1)*Path.biftest_value(end)<=0
                bif_found=1;
                bif_type = 2; % 2: additional; 1: fold bif.; 0: branch point bif; NaN: unknown
                Bifurcation.bif = [Bifurcation.bif,[numel(Path.l_all);bif_type]];
                Bifurcation.flag = 1;
                Bifurcation.scaling = [Bifurcation.scaling,1];
            end
        elseif Opt.bifurcation.determine || Opt.bifurcation.trace || Opt.bifurcation.parameter_trace
            %% determine bifurcation-points:
            Opt_bif = Opt;
            Opt_bif.jacobian = false;
            [bif_solver,default_bif_solver_output] = continuation.solver(Opt_bif,0);
            det_solver_jacobian_red = det(last_jacobian_red);
            Bifurcation.scaling = [Bifurcation.scaling,1/det_solver_jacobian_red];
            if Path.biftest_value(end-1)*Path.biftest_value(end)<=0
                bif_found=1;
                residual_bif = @(x) bifurcation.residual_additional_testfunction(func,x,Opt,Jacobian,Path,Info);
                %% find exact point:
                nds = 11;
                dss = linspace(Path.s_all(end-1)-Path.s_all(end),0,nds);
                dss = dss([6,7,5,8,4,9,3,10,2,11,1]);
                ind_bif = length(Path.l_all);
                bif_type = NaN;
                for i=1:nds
                    dsp = dss(i);
                    [var_bif_predictor,l_bif_predictor] = continuation.predictor(Path,dsp,last_jacobian_red,func,res_corr,Solver,Opt);
                    dscale = aux.get_dscale(Opt,struct('var_all',var_bif_predictor,'l_all',l_bif_predictor));
                    [x_bif,fun_bif,bif_solver_exitflag,bif_solver_output,bif_solver_jacobian] = bif_solver(residual_bif,[var_bif_predictor;l_bif_predictor],dscale);
                    if bif_solver_exitflag>0
%                         Path.s_all = [Path.s_all(1:end-1),Path.s_all(end-1)+[norm(x_bif-[Path.var_all(:,end-1);Path.l_all(end-1)]),norm(x_bif-[Path.var_all(:,end-1);Path.l_all(end-1)])+norm([Path.var_all(:,end);Path.l_all(end)]-x_bif)]];
%                         Path.var_all = [Path.var_all(:,1:end-1),x_bif(1:end-1),Path.var_all(:,end)];
%                         Path.l_all = [Path.l_all(1:end-1),x_bif(end),Path.l_all(end)];
%                         Path.biftest_value=[Path.biftest_value(1:end-1) Opt.bif_additional_testfunction(func,x_bif,Jacobian,Path,Info),Path.biftest_value(end)];
                        % 
                        % get type of bifurcation
                        bif_type = 2; % 2: additional; 1: fold bif.; 0: branch point bif; NaN: unknown
                        break;
                    end
                end
                Bifurcation.bif = [Bifurcation.bif,[ind_bif;bif_type]];
                Bifurcation.flag = 1;
            end
        else
            %% error
            error('No such bifurcation-mode');
        end
    end
end