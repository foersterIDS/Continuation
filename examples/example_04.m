clear;
close all;
clc;
%% residual function:
fun = @(var,l) var-l*sin(10/(l-2.5));
% It is possible that fun returns R or [R,J] with R being the residual and
% J the Jacobian with respect to var or [var;l].
%% settings:
var0 = 0; % initial value of var
l_start = 0; % initial value of l
l_end = 2; % end value of l
ds0 = 0.01; % initial step size
ds_min = 1*10^-5;
ds_max = 2*10^-1;
ssc = {'angle_change';
       'angle_custom';
       'contraction';
       'error';
       'error_alt';
       'fayezioghani';
       'fix';
       'iterations_exponential';
       'iterations_polynomial';
       'multiplicative'; % inital value
       'multiplicative_alt';
       'pid_custom';
       'pid_valli';
       'szyszkowski';
       'yoon'};
n_ssc = numel(ssc);
alpha_reverse = 2.5;
%% run continuation:

% basic options & plot on (alpha_reverse pi, predictor tangential, ds_min/max)
[var_all_0,l_all_0,exitflag_0] = continuation(fun,var0,l_start,l_end,ds0,'plot','on',...
                                                                         'alpha_reverse',alpha_reverse,...
                                                                         'predictor','tangential',...
                                                                         'ds_min',ds_min,'ds_max',ds_max);
drawnow;

% loop over all step size control methods
clc;
figure(2);
clf;
ax = axis;
xlim([l_start,l_end]);
ylim([-2,2]);
xlabel('$\lambda$','interpreter','latex');
ylabel('$v$','interpreter','latex');
grid on;
box on;
drawnow;
for ii=1:n_ssc
    fprintf('run with scc: %s',ssc{ii});
    tic;
    [var_all{ii},l_all{ii},exitflag(ii)] = continuation(fun,var0,l_start,l_end,ds0,'plot','off',...
                                                                                   'display','off',...
                                                                                   'alpha_reverse',alpha_reverse,...
                                                                                   'predictor','tangential',...
                                                                                   'ds_min',ds_min,'ds_max',ds_max,...
                                                                                   'step_size_control',ssc{ii});
    t_ssc(ii) = toc;
    fprintf(' (duration: %.3f s, target reached: %d)\n',t_ssc(ii),l_all{ii}(end)>=l_end);
    hold on;
    plot(l_all{ii},var_all{ii},'-','LineWidth',2,'Color',plot.get_RGB(ii,n_ssc,1));
    hold off;
    drawnow;
end

% result:
[t_ssc_min,ii_min] = min(t_ssc);
fprintf('- - - - - - - - - - - - - - - - - - - - - - - - - - - -\nFastest method: %s (%.3f s)\n',ssc{ii_min},t_ssc_min);