%% path continuation - aux.ison
%
%   Institute of Dynamics and Vibration Research
%   Leibniz University Hannover
%   26.05.2020 - Alwin Förster
%
function [state] = ison(OptSubStruct)
    % Returns true/false weather a sub-struct has one option 'on'
    if isstruct(OptSubStruct)
        state = (sum(structfun(@(x) x,OptSubStruct))>0);
    elseif islogical(OptSubStruct)
        state = OptSubStruct;
    else
        error('Input must be struct or logical!');
    end
end