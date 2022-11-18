%% path continuation - aux.getDscale
%
%   Institute of Dynamics and Vibration Research
%   Leibniz University Hannover
%   06.11.2020 - Alwin Förster
%
function [dscale] = getDscale(Opt,Path)
    if Opt.scaling.dynamicdscale
        % get latest solution
        xLatest = [Path.varAll(:,end); Path.lAll(end)];
        
        % find new dscale
        if numel(Opt.dscaleMin) == 1 || numel(Opt.dscaleMin) == numel(Path.varAll(:,end))+1
            minimal = Opt.dscaleMin;
        else
            error('Size of dscale must either be 1 or equal to size of var0+1!');
        end
        
        dscale = max(abs(xLatest), minimal);
    elseif Opt.scaling.staticdscale
        dscale = Opt.dscale0;
    else
        dscale = ones(size(Opt.dscale0));
    end
end