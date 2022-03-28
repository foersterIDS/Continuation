%% path continuation - homotopy
%
%   Institute of Dynamics and Vibration Research
%   Leibniz University Hannover
%   14.01.2022 - Alwin Förster
%
%% This file is part of continuation.
% 
% If you use continuation, please refer to:
%   A. Förster, foerster@ids.uni-hannover.de
% 
% COPYRIGHT AND LICENSING: 
% continuation Copyright (C) 2022  Alwin Förster
%   							   (foerster@ids.uni-hannover.de)
%                     			   Leibnitz University Hannover
% This program comes with ABSOLUTELY NO WARRANTY. 
% NLvib is free software, you can redistribute and/or modify it under the
% GNU General Public License as published by the Free Software Foundation,
% either version 3 of the License, or (at your option) any later version.
% For details on license and warranty, see http://www.gnu.org/licenses
% or gpl-3.0.txt.
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