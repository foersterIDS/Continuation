%% path continuation - bifurcation.trace
%
%   Institute of Dynamics and Vibration Research
%   Leibniz University Hannover
%   28.10.2020 - Alwin Förster
%   02.07.2021 - Tido Kubatschek
%
function [Path,Bifurcation] = trace(Opt,Opt_is_set,Path,Bifurcation,Solver,Info,func,res_corr)
    bif_trace = Bifurcation.bif(:,Bifurcation.bif(2,:) == 0);
    nbif = numel(bif_trace(1,:));
    Opt_sphere = Opt;
    Opt_sphere = aux.seton(Opt_sphere,'corrector','sphere');
    Opt_trace = Opt;
    Opt_is_set_trace = Opt_is_set;
    Opt_trace.stop_on_crossing = true;
    Opt_is_set_trace.stop_on_crossing = true;
    Opt_trace = aux.seton(Opt_trace,'bifurcation','mark');
    Opt_is_set_trace.bifurcation = true;
    for i=1:nbif
        ind_bif = bif_trace(1,i);
        xdirs_old = [];
        xdirs_trace = [];
        ds_bif = mean(diff(Path.s_all(bif_trace(1,i)+(-1:1))));
        x0 = [Path.var_all(:,bif_trace(1,i));Path.l_all(bif_trace(1,i))];
        residual_bif_sphere = @(x) aux.merge_residuals(func,continuation.corrector(func,Opt_sphere),x,x0,ds_bif,[],Opt_sphere);
        %% find directions of known path
        for j=1:2
            Path_trace = Path;
            if j==1
                Path_trace.var_all = Path.var_all(:,1:bif_trace(1,i));
                Path_trace.l_all = Path.l_all(1:bif_trace(1,i));
                Path_trace.s_all = Path.s_all(1:bif_trace(1,i));
            else
                Path_trace.var_all = Path.var_all(:,end:-1:bif_trace(1,i));
                Path_trace.l_all = Path.l_all(end:-1:bif_trace(1,i));
                Path_trace.s_all = abs(Path.s_all(end:-1:bif_trace(1,i))-Path.s_all(end));
            end
            [var_bif_predictor,l_bif_predictor] = continuation.predictor(Path_trace,ds_bif,[],func,res_corr,Solver,Opt_sphere);
            x_bif_predictor = [var_bif_predictor;l_bif_predictor];
            dscale = aux.get_dscale(Opt,struct('var_all',var_bif_predictor,'l_all',l_bif_predictor));
            [x_bif_ij,~,solver_bif_exitflag] = Solver.main(residual_bif_sphere,x_bif_predictor,dscale);
            if solver_bif_exitflag>0
                xdirs_old = [xdirs_old,x_bif_ij-x0];
                residual_bif_sphere = @(x) aux.deflation(residual_bif_sphere,x_bif_ij,x,Opt_sphere);
            end
        end
        %% unknown paths
        if ~isempty(xdirs_old)
            %% find directions of unknown paths by using random direction
            if Opt.bif_rand_dir
                for j=1:Opt.n_bif_search
                    dx_bif_predictor = randn(numel(x0),1);
                    for ki = 1:2
                        x_bif_predictor = x0+(-1)^ki*ds_bif*dx_bif_predictor/norm(dx_bif_predictor);
                        dscale = aux.get_dscale(Opt,struct('var_all',x_bif_predictor(1:end-1,:),'l_all',x_bif_predictor(end,:)));
                        [x_bif_ij,~,solver_bif_exitflag] = Solver.main(residual_bif_sphere,x_bif_predictor,dscale);
                        if solver_bif_exitflag>0 && norm(x_bif_ij-x0)>=ds_bif*0.99 && norm(x_bif_ij-x0)<=ds_bif*1.01
                            xdirs_trace = [xdirs_trace,x_bif_ij-x0];
                            residual_bif_sphere = @(x) aux.deflation(residual_bif_sphere,x_bif_ij,x,Opt_sphere);
                        end
                    end              
                end
            end
            %% use tangent vectors - experimental!
            for j=1:(numel(Bifurcation.dirs)/2)
                [bif_num, bif_dir] = Bifurcation.dirs{j,:};
                if bif_num == ind_bif
                    if ~isempty(bif_dir)
                        for ki = 1:2
                            x_bif_predictor = x0+(-1)^ki*ds_bif*bif_dir/norm(bif_dir);
                            dscale = aux.get_dscale(Opt,struct('var_all',x_bif_predictor(1:end-1,:),'l_all',x_bif_predictor(end,:)));
                            [x_bif_ij,~,solver_bif_exitflag] = Solver.main(residual_bif_sphere,x_bif_predictor,dscale);
                            if solver_bif_exitflag>0 && norm(x_bif_ij-x0)>=ds_bif*0.99 && norm(x_bif_ij-x0)<=ds_bif*1.01
                                xdirs_trace = [xdirs_trace,x_bif_ij-x0];
                                residual_bif_sphere = @(x) aux.deflation(residual_bif_sphere,x_bif_ij,x,Opt_sphere);
                            end
                        end
                    end
                end
            end
            %% trace unknown paths
            ntrace = numel(xdirs_trace)/numel(x0);
            for j=1:ntrace
                Opt_trace.direction = xdirs_trace(:,j)/norm(xdirs_trace(:,j));
                Opt_is_set_trace.direction = true;
                Opt_trace.l_0 = x0(end);
                Opt_is_set_trace.l_0 = true;
                Opt_trace.initial_deflation_points = x0+xdirs_old;
                Opt_is_set_trace.initial_deflation_points = true;
                Opt_trace.Opt_is_set = Opt_is_set_trace; % (set Opt_is_set with Opt)
                [var_j,l_j,exitflag_j,Bifurcation_j,s_j] = continuation(func,x0(1:end-1),Info.l_start,Info.l_end,ds_bif,'opt',Opt_trace);
                if exitflag_j>0
                    Path.var_all = [Path.var_all,NaN(numel(x0)-1,1),var_j];
                    Path.l_all = [Path.l_all,NaN,l_j];
                    Path.s_all = [Path.s_all,NaN,s_j];
                    Bifurcation.bif = [Bifurcation.bif,Bifurcation_j.bif];
                end
            end
        end
    end
end