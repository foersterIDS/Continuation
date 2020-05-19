%% path continuation - arclength
%
%   Institute of Dynamics and Vibration Research
%   Leibniz University Hannover
%   08.05.2020 - Alwin Förster
%
function [dsn] = arclength(ds,ds0,error_counter)
    if error_counter==0
        dsn = min([ds0,ds*2]);
    else
        dsn = ds/2;
    end
end