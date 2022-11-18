%% path continuation - aux.descaleValues
%
%   Institute of Dynamics and Vibration Research
%   Leibniz University Hannover
%   03.11.2020 - Tido Kubatschek
%
function [xDescaled] = descaleValues(Opt, varAll, lAll)
    if ison(Opt.scaling)
        xAll = [varAll; lAll];
        Descale = diag(Opt.dscale);
        xDescaled = Descale * xAll(:,end);
    end
end

