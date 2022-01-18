%% path continuation - aux.setoff
%
%   Institute of Dynamics and Vibration Research
%   Leibniz University Hannover
%   18.01.2022 - Alwin Förster
%
function [Opt] = setoff(Opt,substruct)
    if isfield(Opt,substruct)
        fields = fieldnames(Opt.(substruct));
        for i=1:numel(fields)
            Opt.(substruct).(fields{i}) = false;
        end
    else
        error('no such substruct');
    end
end