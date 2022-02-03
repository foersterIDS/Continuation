%% path continuation - aux.clear_struct
%
%   Institute of Dynamics and Vibration Research
%   Leibniz University Hannover
%   03.02.2022 - Alwin Förster
%
function [struct_out] = clear_struct(struct_in)
    fn = fieldnames(struct_in);
    nfields = numel(fn);
    struct_out = cell2struct(cell(nfields,1),fn(:),1);
    for ii=1:nfields
        if isstruct(struct_in.(fn{ii}))
            struct_out.(fn{ii}) = aux.clear_struct(struct_in.(fn{ii}));
        end
    end
end