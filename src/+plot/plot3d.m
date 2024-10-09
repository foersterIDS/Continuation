function [x,y,z] = plot3d(l,var,NameValueArgs)
    arguments
        l (1,1) double
        var (:,1) double
        NameValueArgs.idxY (1,1) double {mustBeInteger(NameValueArgs.idxY),mustBeGreaterThan(NameValueArgs.idxY,0)}
        NameValueArgs.idxZ (1,1) double {mustBeInteger(NameValueArgs.idxZ),mustBeGreaterThan(NameValueArgs.idxZ,0)}
    end
    x = l;
    nVar = length(var);
    if nVar < 2
        error('Cannot use plot3d if nVar < 2')
    elseif nVar > 2
        if ~(isempty(NameValueArgs.idxY) && isempty(NameValueArgs.idxZ))
            if any([NameValueArgs.idxY,NameValueArgs.idxZ] > nVar)
                error('idxY and idxZ cannot be greater than nVar.');
            end
            y = var(NameValueArgs.idxY,end);
            z = var(NameValueArgs.idxZ,end);
        else
            warning('Using first two entries of varAll');
            y = var(1,end);
            z = var(2,end);
        end
    else
        y = var(1,end);
        z = var(2,end);
    end
end