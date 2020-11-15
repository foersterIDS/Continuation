%% path continuation - get_dscale
%
%   Institute of Dynamics and Vibration Research
%   Leibniz University Hannover
%   06.11.2020 - Alwin Förster
%
function [dscale] = get_dscale(Opt,var_all,l_all)
    if Opt.scaling.dynamicdscale
        % get latest solution
        x_latest = [var_all(:,end); l_all(end)];
        
        % save latest dscale
        dscaleold = Opt.dscale;
        
        % find new dscale
        Opt.dscale(1:end-1) = max(abs(x_latest(1:end-1)),...
            Opt.dscale0(1:end-1));
        
        dscale = Opt.dscale;
        
    elseif Opt.scaling.staticdscale
        dscale = Opt.dscale0;
    else
        dscale = ones(size(Opt.dscale0));
    end
end