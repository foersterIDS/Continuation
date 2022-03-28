clear;
close all;
clc;
%% residual function:
fun = @(var,l) var-l-sin(3/2*pi*l);
% It is possible that fun returns R or [R,J] with R being the residual and
% J the Jacobian with respect to var or [var;l].
%% settings:
var0 = 0; % initial value of var
l_start = 0; % initial value of l
l_end = 2; % end value of l
ds0 = 0.01; % initial step size
%% run continuation:

% basic options
[var_all,l_all,exitflag] = continuation(fun,var0,l_start,l_end,ds0);
input('Press any key to continue...');

% use plot and max. step size
[var_all,l_all,exitflag] = continuation(fun,var0,l_start,l_end,ds0,...
                                        'plot','on',... % show live plot
                                        'ds_max',0.25); % max. step size
input('Press any key to continue...');

% use deteiled plot and max. step size
[var_all,l_all,exitflag] = continuation(fun,var0,l_start,l_end,ds0,...
                                        'plot','detail',... % show detailed live plot
                                        'ds_max',0.25); % max. step size
input('Press any key to continue...');

% use plot, max. step size and corrector 'orthogonal' instead of 'sphere'
[var_all,l_all,exitflag] = continuation(fun,var0,l_start,l_end,ds0,...
                                        'plot','on',... % show live plot
                                        'ds_max',0.25,... % max. step size
                                        'corrector','orthogonal'); % corrector orthogognal
input('Press any key to continue...');

% use plot, max. step size and polynomial predictor of order 3 instead of a
% tangential one
[var_all,l_all,exitflag] = continuation(fun,var0,l_start,l_end,ds0,...
                                        'plot','on',... % show live plot
                                        'ds_max',0.25,... % max. step size
                                        'predictor','polynomial',... % polynomial predictor
                                        'predictor_polynomial_order',3); % polynomial order