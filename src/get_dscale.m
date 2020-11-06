%% path continuation - get_dscale
%
%   Institute of Dynamics and Vibration Research
%   Leibniz University Hannover
%   06.11.2020 - Alwin Förster
%
function [dscale] = get_dscale(Opt,var_all,l_all)
    if Opt.scaling.dynamicdscale
        % TODO
        dscale = Opt.dscale0;
    elseif Opt.scaling.staticdscale
        dscale = Opt.dscale0;
    else
        dscale = Opt.dscale0;
    end
end