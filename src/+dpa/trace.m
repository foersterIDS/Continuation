%% path continuation - dpa.trace
%
%   Institute of Dynamics and Vibration Research
%   Leibniz University Hannover
%   04.04.2022 - Alwin Förster
%
function [Path] = trace(fun,dpa_points,Info,Opt,Path)
    %% settings for dpa:
    dpa_trace = dpa_points;
    npoints = numel(dpa_trace);
    Opt_trace = Opt;
    Opt_trace.direction = [zeros(size(Opt_trace.direction));sign(Opt.g_target-Opt.g_0)];
    Opt_trace.l_0 = Opt.g_0;
    Opt_trace.l_target = Opt.g_target;
    Opt_trace.dpa_gamma_var = true;
    Path_dpa = struct('var_all',[],'l_all',[],'s_all',[]);
    %% start dpa:
    for ii=1:npoints
        i = dpa_trace(ii);
        func = @(x,g) dpa.merge_residuals(Opt,fun,Opt.dpa_residual,x,g);
        x0 = [Path.var_all(:,i);Path.l_all(i)];
        ds_dpa = max([mean(diff(Path.s_all(i+(-1:1))))/100,Opt.ds_min]);
        Opt_trace.dscale0 = max(abs([x0;Opt.g_0]),10^-8*ones(numel(x0)+1,1));
        [var_i,l_i,~,~,s_i] = continuation(func,x0,Opt.g_0,Opt.g_target,ds_dpa,'Opt',Opt_trace);
        Path_dpa.var_all = [Path_dpa.var_all,NaN(Info.nv+1,1),var_i];
        Path_dpa.l_all = [Path_dpa.l_all,NaN,l_i];
        Path_dpa.s_all = [Path_dpa.s_all,NaN,s_i];
    end
    %% output:
    Path_temp = Path;
    if ~isempty(Path_dpa.var_all)
        Path.var_all = [Path_temp.var_all,Path_dpa.var_all(1:Info.nv,:)];
        Path.l_all = [Path_temp.l_all,Path_dpa.var_all(Info.nv+1,:);
                      Opt.g_0*ones(size(Path_temp.l_all)),Path_dpa.l_all];
        Path.s_all = [Path_temp.s_all,Path_dpa.s_all];
    else
        Path.l_all = [Path_temp.l_all;Opt.g_0*ones(size(Path_temp.l_all))];
    end
end