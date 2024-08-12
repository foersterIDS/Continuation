%% path continuation - aux.getDscale
%
%   Institute of Dynamics and Vibration Research
%   Leibniz University Hannover
%   06.11.2020 - Alwin Förster
%
function [dscale] = getDscale(oih,xScale)
    arguments
        oih (1,1) aux.OptInfoHandle
        xScale (:,1) double = oih.path.xAll(:,end)
    end
    if oih.opt.scaling.dynamicdscale
        % find new dscale
        if numel(oih.opt.dscaleMin) == 1 || numel(oih.opt.dscaleMin) == numel(oih.path.varAll(:,end))+1
            minimal = oih.opt.dscaleMin;
        else
            error('Size of dscale must either be 1 or equal to size of var0+1!');
        end
        
        dscale = max(abs(xScale), minimal);
    elseif oih.opt.scaling.staticdscale
        dscale = oih.opt.dscale0;
    else
        dscale = ones(size(oih.opt.dscale0));
    end
end