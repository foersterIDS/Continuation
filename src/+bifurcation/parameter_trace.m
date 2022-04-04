%% path continuation - bifurcation.parameter_trace
%
%   Institute of Dynamics and Vibration Research
%   Leibniz University Hannover
%   31.03.2022 - Alwin Förster
%
function [Path] = parameter_trace(Opt,Path,Bifurcation,Info,fun)
    %% settings for dpa:
    bif_trace = Bifurcation.bif(1,:);
    nbif = numel(bif_trace);
    Opt_trace = Opt;
    Opt_trace = aux.setoff(Opt_trace,'bifurcation');
    Opt_trace.direction = [zeros(size(Opt_trace.direction));sign(Opt.g_target-Opt.g_0)];
    Opt_trace.l_0 = Opt.g_0;
    Opt_trace.l_target = Opt.g_target;
    Opt_trace.dpa_gamma_var = true;
    Path_bifs = struct('var_all',[],'l_all',[],'s_all',[]);
    %% start dpa:
    for ii=1:nbif
        i = bif_trace(ii);
        sc = Bifurcation.scaling(ii);
        res_dpa = @(v,l,g) dpa.res_bif(fun,[v;l],g,Opt,sc);
        func = @(x,g) dpa.merge_residuals(Opt,fun,res_dpa,x,g);
        x0 = [Path.var_all(:,i);Path.l_all(i)];
        ds_bif = mean(diff(Path.s_all(i+(-1:1))));
        Opt_trace.dscale0 = max(abs([x0;Opt.g_0]),10^-8*ones(numel(x0)+1,1));
        [var_i,l_i,~,~,s_i] = continuation(func,x0,Opt.g_0,Opt.g_target,ds_bif,'Opt',Opt_trace);
        Path_bifs.var_all = [Path_bifs.var_all,NaN(Info.nv+1,1),var_i];
        Path_bifs.l_all = [Path_bifs.l_all,NaN,l_i];
        Path_bifs.s_all = [Path_bifs.s_all,NaN,s_i];
    end
    %% output:
    Path_temp = Path;
    if ~isempty(Path_bifs.var_all)
        Path.var_all = [Path_temp.var_all,Path_bifs.var_all(1:Info.nv,:)];
        Path.l_all = [Path_temp.l_all,Path_bifs.var_all(Info.nv+1,:);
                      Opt.g_0*ones(size(Path_temp.l_all)),Path_bifs.l_all];
        Path.s_all = [Path_temp.s_all,Path_bifs.s_all];
    else
        Path.l_all = [Path_temp.l_all;Opt.g_0*ones(size(Path_temp.l_all))];
    end
end