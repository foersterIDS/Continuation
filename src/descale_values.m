%% path continuation - rescale_values
%
%   Institute of Dynamics and Vibration Research
%   Leibniz University Hannover
%   03.11.2020 - Tido Kubatschek
%
function [x_descaled] = descale_values(Opt, var_all, l_all)
    if ison(Opt.scaling)
        x_all = [var_all; l_all];
        Descale = diag(Opt.Dscale);
        x_descaled = Descale * x_all(:,end);
    end
end

