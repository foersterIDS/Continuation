%% path continuation - homotopy
%
%   Institute of Dynamics and Vibration Research
%   Leibniz University Hannover
%   14.01.2022 - Alwin Förster
%
function [ v_solution, exitflag, var_all, l_all ] = homotopy( fun, var0, varargin )
    %% initialize
    %
    warning on;
    Opt = continuation.input(varargin,fun,var0);
    %
    %% homotopy
    %
    [v_solution, exitflag, var_all, l_all] = homotopy.h_continuation(fun,var0,Opt);
    %
end