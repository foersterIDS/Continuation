%% path continuation - trace_bifurcations
%
%   Institute of Dynamics and Vibration Research
%   Leibniz University Hannover
%   28.10.2020 - Alwin Förster
%   02.07.2021 - Tido Kubatschek
%
function [var_all,l_all,s_all,bif] = trace_bifurcations(Opt,var_all,l_all,s_all,bif,solver,fun,l_start,l_end,res_corr,predictor_solver,bif_dirs)
    nbif = numel(bif(2,bif(2,:) == 0));
    bif_trace = bif(:,bif(2,:) == 0);
    Opt_sphere = Opt;
    Opt_sphere = seton(Opt_sphere,'corrector','sphere');
    Opt_trace = Opt;
    Opt_trace.stop_on_bifurcation = true;
    Opt_trace = seton(Opt_trace,'bifurcation','determine');
    for i=1:nbif
        ind_bif = bif_trace(1,i);
        xdirs_old = [];
        xdirs_trace = [];
        ds_bif = mean(diff(s_all(bif_trace(1,i)+(-1:1))));
        x0 = [var_all(:,bif_trace(1,i));l_all(bif_trace(1,i))];
        residual_bif_sphere = @(x) merge_residuals(fun,residual_corrector(Opt_sphere),x,x0,ds_bif,Opt_sphere);
        %% find directions of known path
        for j=1:2
            [var_bif_predictor,l_bif_predictor] = predictor(var_all(:,1:bif_trace(1,i)),l_all(1:bif_trace(1,i)),s_all(1:bif_trace(1,i)),(-1)^j*ds_bif,[],fun,res_corr,predictor_solver,Opt_sphere);
            x_bif_predictor = [var_bif_predictor;l_bif_predictor];
            dscale = get_dscale(Opt,var_bif_predictor,l_bif_predictor);
            [x_bif_ij,~,solver_bif_exitflag] = solver(residual_bif_sphere,x_bif_predictor,dscale);
            if solver_bif_exitflag>0
                xdirs_old = [xdirs_old,x_bif_ij-x0];
                residual_bif_sphere = @(x) deflation(residual_bif_sphere,x_bif_ij,x,Opt_sphere);
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
                        dscale = get_dscale(Opt,x_bif_predictor(1:end-1,:),x_bif_predictor(end,:));
                        [x_bif_ij,~,solver_bif_exitflag] = solver(residual_bif_sphere,x_bif_predictor,dscale);
                        if solver_bif_exitflag>0 && norm(x_bif_ij-x0)>=ds_bif*0.99 && norm(x_bif_ij-x0)<=ds_bif*1.01
                            xdirs_trace = [xdirs_trace,x_bif_ij-x0];
                            residual_bif_sphere = @(x) deflation(residual_bif_sphere,x_bif_ij,x,Opt_sphere);
                        end
                    end              
                end
            end
            %% use tangent vectors - experimental!
            for j=1:(numel(bif_dirs)/2)
                [bif_num, bif_dir] = bif_dirs{j,:};
                if bif_num == ind_bif
                    if ~isempty(bif_dir)
                        for ki = 1:2
                            x_bif_predictor = x0+(-1)^ki*ds_bif*bif_dir/norm(bif_dir);
                            dscale = get_dscale(Opt,x_bif_predictor(1:end-1,:),x_bif_predictor(end,:));
                            [x_bif_ij,~,solver_bif_exitflag] = solver(residual_bif_sphere,x_bif_predictor,dscale);
                            if solver_bif_exitflag>0 && norm(x_bif_ij-x0)>=ds_bif*0.99 && norm(x_bif_ij-x0)<=ds_bif*1.01
                                xdirs_trace = [xdirs_trace,x_bif_ij-x0];
                                residual_bif_sphere = @(x) deflation(residual_bif_sphere,x_bif_ij,x,Opt_sphere);
                            end
                        end
                    end
                end
            end
            %% trace unknown paths
            ntrace = numel(xdirs_trace)/numel(x0);
            for j=1:ntrace
                Opt_trace.direction = xdirs_trace(:,j)/norm(xdirs_trace(:,j));
                Opt_trace.l_0 = x0(end);
                [var_j,l_j,exitflag_j,bif_j,s_j] = continuation(fun,x0(1:end-1),l_start,l_end,ds_bif,'opt',Opt_trace);
                if exitflag_j>0
                    var_all = [var_all,NaN(numel(x0)-1,1),var_j];
                    l_all = [l_all,NaN,l_j];
                    s_all = [s_all,NaN,s_j];
                    bif = [bif,bif_j];
                end
            end
        end
    end
end