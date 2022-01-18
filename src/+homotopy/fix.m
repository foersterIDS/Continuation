%% path continuation - homotopy.fix
%
%   Institute of Dynamics and Vibration Research
%   Leibniz University Hannover
%   19.05.2020 - Alwin Förster
%
function [f,J] = fix(x,x0)
    f = x-x0;
    J = ones(size(x0));
end