%% path continuation - aux.ison
%
%   Institute of Dynamics and Vibration Research
%   Leibniz University Hannover
%   26.05.2020 - Alwin Förster
%
function [state,type] = ison(optSubStruct)
    % Returns true/false weather a sub-struct has one option 'on'
    if isstruct(optSubStruct)
        onTemp = structfun(@(x) x,optSubStruct);
        state = sum(onTemp)>0;
        fn = fieldnames(optSubStruct);
        if state
            type = fn{onTemp};
        else
            type = [];
        end
    elseif islogical(optSubStruct)
        state = optSubStruct;
        type = [];
    else
        error('Input must be struct or logical!');
    end
end