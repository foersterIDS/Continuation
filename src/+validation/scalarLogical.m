%% path continuation - validation.scalarLogical
%
%   Institute of Dynamics and Vibration Research
%   Leibniz University Hannover
%   26.05.2023 - Alwin Förster
%
function [valid] = scalarLogical(input)
    if ischar(input) || isstring(input)
        if strcmpi(input,'on')
            valid = true;
        elseif strcmpi(input,'off')
            valid = true;
        else
            valid = false;
        end
    elseif islogical(logical(input)) && (numel(logical(input))==1)
        valid = true;
    else
        valid = false;
    end
end