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
input('Press any key to continue...');

%% dpa residual function:
fun_dpa = @(var,l,g) (l-2*var-g)*(l-(var-2*g).^2-5);
% It is possible that fun returns R or [R,J] with R being the residual and
% J the Jacobian with respect to var or [var;l].
%% dpa settings:
l_end = 25; % end value of l
g_0 = 1; % initial value of g
g_target = 3; % end value of g
ds0 = 0.001; % initial step size
ds_max = 0.05; % max. step size
%% run continuation:

% show dpa plot, use ds_max and set bifurcation detection to parameter_trace
[var_all_dpa,l_all_dpa,exitflag_dpa] = continuation(fun_dpa,var0,l_start,l_end,ds0,'plot','dpa',...
                                                                                   'ds_max',ds_max,...
                                                                                   'bifurcation','parameter_trace',...
                                                                                   'g_0',g_0,...
                                                                                   'g_target',g_target);
input('Press any key to continue...');

% show plot, use ds_max and set bifurcation detection to trace, fix at g_0
[var_all_trace_g_0,l_all_trace_g_0,exitflag_trace_g_0] = continuation(@(v,l) fun_dpa(v,l,g_0),var0,l_start,l_end,ds0,'plot','on',...
                                                                                                                     'ds_max',ds_max,...
                                                                                                                     'bifurcation','trace');
input('Press any key to continue...');

% show plot, use ds_max and set bifurcation detection to trace, fix at g_target
[var_all_trace_g_target,l_all_trace_g_target,exitflag_trace_g_target] = continuation(@(v,l) fun_dpa(v,l,g_target),var0,l_start,l_end,ds0,'plot','on',...
                                                                                                                                         'ds_max',ds_max,...
                                                                                                                                         'bifurcation','trace');
input('Press any key to continue...');

%% plot:
figure(8);
clf;
plot3(l_all_dpa(1,:),l_all_dpa(2,:),var_all_dpa,'k-','LineWidth',2);
hold on;
plot3(l_all_trace_g_0,g_0*ones(size(l_all_trace_g_0)),var_all_trace_g_0,'b-','LineWidth',2);
plot3(l_all_trace_g_target,g_target*ones(size(l_all_trace_g_target)),var_all_trace_g_target,'r-','LineWidth',2);
hold off;
grid on;
xlim([l_start,l_end]);
ylim([g_0,g_target]);
drawnow;