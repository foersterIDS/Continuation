%% path continuation - seton
%
%   Institute of Dynamics and Vibration Research
%   Leibniz University Hannover
%   28.10.2020 - Alwin Förster
%
function [Opt] = seton(Opt,substruct,option)
    if isfield(Opt,substruct)
        if isfield(Opt.(substruct),option)
            fields = fieldnames(Opt.(substruct));
            for i=1:numel(fields)
                Opt.(substruct).(fields{i}) = false;
            end
            Opt.(substruct).(option) = true;
        else
            error('no such option');
        end
    else
        error('no such substruct');
    end
end