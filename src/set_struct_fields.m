%% path continuation - set_struct_fields
%
%   Institute of Dynamics and Vibration Research
%   Leibniz University Hannover
%   10.01.2022 - Tido Kubatschek
%
function [ s_struct ] = set_struct_fields(s_struct, varargin)
    var_len = length(varargin);
    
    if mod(var_len,2) ~= 0
        error('There must be an even number of fieldnames and values!\n')
    end
    
    for k = 1:2:(var_len-1)
        field_name = varargin{k};
        
        if ~ischar(field_name)
            error('Input is not a string!');
        end
        
        field_value = varargin{k+1};
        s_struct = setfield(s_struct,field_name, field_value);
    end
    
end

