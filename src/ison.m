%% path continuation - ison
%
%   Institute of Dynamics and Vibration Research
%   Leibniz University Hannover
%   26.05.2020 - Alwin Förster
%
function [state] = ison(Opt_sub_struct)
    % Returns true/false weather a sub-struct has one option 'on'
    if isstruct(Opt_sub_struct)
        state = (sum(structfun(@(x) x,Opt_sub_struct))>0);
    elseif islogical(Opt_sub_struct)
        state = Opt_sub_struct;
    else
        error('Input must be struct or logical!');
    end
end