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
        
        % find new dscale
        if numel(Opt.dscale_min) == 1 || numel(Opt.dscale_min) == numel(var_all(:,end))+1
            minimal = Opt.dscale_min;
        else
            error('Size of dscale must either be 1 or equal to size of var0+1!');
        end
        
        dscale = max(abs(x_latest), minimal);
    elseif Opt.scaling.staticdscale
        dscale = Opt.dscale0;
    else
        dscale = ones(size(Opt.dscale0));
    end
end