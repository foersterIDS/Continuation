clear;
close all;
clc;
%% residual function:
fun = @(var,l) var.^2+l.^2-1;
% It is possible that fun returns R or [R,J] with R being the residual and
% J the Jacobian with respect to var or [var;l].
%% settings:
var0 = 1/sqrt(2); % initial value of var
l_start = -2; % left boundary of l
l_end = +2; % right boundary of l
ds0 = 0.001; % initial step size
l_0 = 1/sqrt(2); % initial value of l
%% run continuation:

% start at l_0 and show plot
[var_all,l_all,exitflag] = continuation(fun,var0,l_start,l_end,ds0,'l_0',l_0,...
                                                                   'plot','on');
input('Press any key to continue...');

% start at l_0 and show plot, turn closed_curve_detection off and
% set n_step_max to 500
[var_all,l_all,exitflag] = continuation(fun,var0,l_start,l_end,ds0,'l_0',l_0,...
                                                                   'plot','on',...
                                                                   'closed_curve_detection','off',...
                                                                   'n_step_max',500);