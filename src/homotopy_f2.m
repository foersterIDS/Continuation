%% path continuation - homotopy_f2
%
%   Institute of Dynamics and Vibration Research
%   Leibniz University Hannover
%   04.06.2020 - Alwin Förster
%
function [f] = homotopy_f2(R,x,l)
    f = 1+R(x).^2-l;
end