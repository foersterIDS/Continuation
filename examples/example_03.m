clear;
close all;
clc;
%% residual function:
fun = @(var,l) (l-2*var-1)*(l-(var-2).^2-5);
% It is possible that fun returns R or [R,J] with R being the residual and
% J the Jacobian with respect to var or [var;l].
%% settings:
var0 = 1; % initial value of var
l_start = 0; % initial value of l
l_end = 10; % end value of l
ds0 = 0.01; % initial step size
ds_max = 0.1; % max. step size
%% run continuation:

% show plot and use ds_max
[var_all,l_all,exitflag] = continuation(fun,var0,l_start,l_end,ds0,'plot','on',...
                                                                   'ds_max',ds_max);
input('Press any key to continue...');

% show plot, use ds_max and turn bifurcation detection off
[var_all,l_all,exitflag] = continuation(fun,var0,l_start,l_end,ds0,'plot','on',...
                                                                   'ds_max',ds_max,...
                                                                   'bifurcation','off');
input('Press any key to continue...');

% show plot, use ds_max and set bifurcation detection to mark
[var_all,l_all,exitflag] = continuation(fun,var0,l_start,l_end,ds0,'plot','on',...
                                                                   'ds_max',ds_max,...
                                                                   'bifurcation','mark');
input('Press any key to continue...');

% show plot, use ds_max and set bifurcation detection to determine
[var_all,l_all,exitflag] = continuation(fun,var0,l_start,l_end,ds0,'plot','on',...
                                                                   'ds_max',ds_max,...
                                                                   'bifurcation','determine');
input('Press any key to continue...');

% show plot, use ds_max and set bifurcation detection to trace
[var_all,l_all,exitflag] = continuation(fun,var0,l_start,l_end,ds0,'plot','on',...
                                                                   'ds_max',ds_max,...
                                                                   'bifurcation','trace');