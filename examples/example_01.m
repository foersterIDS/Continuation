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
t0 = tic;
[var_all0,l_all0,exitflag0] = continuation(fun,var0,l_start,l_end,ds0);
t0 = toc(t0);
input('Press any key to continue...');

% use plot and max. step size
t1 = tic;
[var_all1,l_all1,exitflag1] = continuation(fun,var0,l_start,l_end,ds0,...
                                           'plot','on',... % show live plot
                                           'ds_max',0.25); % max. step size
t1 = toc(t1);
input('Press any key to continue...');
close all;

% use deteiled plot and max. step size
t2 = tic;
[var_all2,l_all2,exitflag2] = continuation(fun,var0,l_start,l_end,ds0,...
                                           'plot','detail',... % show detailed live plot
                                           'ds_max',0.25); % max. step size
t2 = toc(t2);
input('Press any key to continue...');
close all;

% use plot, max. step size and corrector 'orthogonal' instead of 'sphere'
t3 = tic;
[var_all3,l_all3,exitflag3] = continuation(fun,var0,l_start,l_end,ds0,...
                                           'plot','on',... % show live plot
                                           'ds_max',0.25,... % max. step size
                                           'corrector','orthogonal'); % corrector orthogognal
t3 = toc(t3);
input('Press any key to continue...');
close all;

% use plot, max. step size and polynomial predictor of order 3 instead of a
% tangential one
t4 = tic;
[var_all4,l_all4,exitflag4] = continuation(fun,var0,l_start,l_end,ds0,...
                                           'plot','on',... % show live plot
                                           'ds_max',0.25,... % max. step size
                                           'corrector','orthogonal',... % corrector orthogognal
                                           'predictor','polynomial',... % polynomial predictor
                                           'predictor_polynomial_degree',3); % polynomial order
t4 = toc(t4);
input('Press any key to continue...');
close all;

%% eval.:
clc;
fprintf('- - - - - - - - - - - - - - - - - - - - - - - - -\n\n');
fprintf('settings (0): %d steps in %.3f s (t/s: %.1f ms)\n',numel(l_all0),t0,10^3*t0/numel(l_all0));
fprintf('settings (1): %d steps in %.3f s (t/s: %.1f ms)\n',numel(l_all1),t1,10^3*t1/numel(l_all1));
fprintf('settings (2): %d steps in %.3f s (t/s: %.1f ms)\n',numel(l_all2),t2,10^3*t2/numel(l_all2));
fprintf('settings (3): %d steps in %.3f s (t/s: %.1f ms)\n',numel(l_all3),t3,10^3*t3/numel(l_all3));
fprintf('settings (4): %d steps in %.3f s (t/s: %.1f ms)\n',numel(l_all4),t4,10^3*t4/numel(l_all4));
fprintf('\n- - - - - - - - - - - - - - - - - - - - - - - - -\n');