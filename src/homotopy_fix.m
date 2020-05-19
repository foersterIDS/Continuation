%% path continuation - homotopy_fix
%
%   Institute of Dynamics and Vibration Research
%   Leibniz University Hannover
%   19.05.2020 - Alwin Förster
%
function [f,J] = homotopy_fix(x,x0)
    f = x-x0;
    J = ones(size(x0));
end