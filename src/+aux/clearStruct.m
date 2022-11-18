%% path continuation - aux.clearStruct
%
%   Institute of Dynamics and Vibration Research
%   Leibniz University Hannover
%   03.02.2022 - Alwin Förster
%
function [structOut] = clearStruct(structIn)
    fn = fieldnames(structIn);
    nfields = numel(fn);
    structOut = cell2struct(cell(nfields,1),fn(:),1);
    for ii=1:nfields
        if isstruct(structIn.(fn{ii}))
            structOut.(fn{ii}) = aux.clearStruct(structIn.(fn{ii}));
        end
    end
end