%% path continuation - homotopy.squared
%
%   Institute of Dynamics and Vibration Research
%   Leibniz University Hannover
%   21.10.2022 - Alwin Förster
%
function [f,J] = squared(x,x0)
    f = (x-x0).^2;
    J = diag(2*(x-x0));
end