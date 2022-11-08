clear;
close all;
clc;
%% residual function:
fun = @(var,l) var-l;
% It is possible that fun returns R or [R,J] with R being the residual and
% J the Jacobian with respect to var or [var;l].
%% settings:
var0 = 0; % initial value of var
l_start = 0; % initial value of l
l_end = 5; % end value of l
ds0 = 0.01; % initial step size
ds_min = 1*10^-5;
ds_max = 2*10^-1;

%% event based stepsize control with two simple events
% create StepSizeControlEvent object
event_obj = StepSizeControlEvent();

% Add first event
% set new ds_max to be 5e-2 when lambda is between 1 and 1.4
event_condition = @(Path) Path.l_all(end) >= 1 && Path.l_all(end) <= 1.4;
event_needed_parameters = {'Path'};
event_ds_min = ds_min;
event_ds_max = 1e-3;
event_counter_max = 1;
event_obj = event_obj.addEvent(event_condition,event_needed_parameters,event_ds_min,event_ds_max,event_counter_max);

% Add second event
% set new ds_max to be 1 when lambda is greater than 3
event_condition = @(Path) Path.l_all(end) > 3;
event_needed_parameters = {'Path'};
event_ds_min = ds_min;
event_ds_max = 1e-2;
event_counter_max = 1;
event_obj = event_obj.addEvent(event_condition,event_needed_parameters,event_ds_min,event_ds_max,event_counter_max);


%% run continuation:

% basic options & plot on (predictor tangential, ds_min/max) & activate
% event based step size control
[var_all,l_all,exitflag,~,s_all] = continuation(fun,var0,l_start,l_end,ds0,'plot','on',...
                                                                         'predictor','tangential',...
                                                                         'ds_min',ds_min,'ds_max',ds_max,...
                                                                         'step_size_event',true,...
                                                                         'event_user_input',event_obj.getEvents);
%% Also plot used stepsize
figure('Units','normalized','Position',[0.2,0.2,0.6,0.6]);

ds = [0,s_all(2:end) - s_all(1:end-1)];
pl(1) = plot(l_all, ds,'b--x',LineWidth=2); hold on;
xlabel('$\lambda$',Interpreter='latex');
ylabel('$\Delta s(\lambda)$',Interpreter='latex');
xl = xlim();
yl = ylim();
pl(2) = fill([1, 1.4, 1.4, 1],[yl(1),yl(1),yl(2),yl(2)],'r','EdgeColor','none','FaceAlpha',0.4); hold on;
pl(3) = fill([3, xl(2), xl(2), 3],[yl(1),yl(1),yl(2),yl(2)],'g','EdgeColor','none','FaceAlpha',0.4);
pl(4) = fill([xl(1), 1, 1, xl(1)],[yl(1),yl(1),yl(2),yl(2)],[0.5,0.5,0.5],'EdgeColor','none','FaceAlpha',0.4);
fill([1.4, 3, 3, 1.4],[yl(1),yl(1),yl(2),yl(2)],[0.5,0.5,0.5],'EdgeColor','none','FaceAlpha',0.4);
uistack(pl(1),'top');
legend(pl(1:4),'used stepsize','event 1', 'event 2', 'no event');
fontsize(gcf, 13, 'points')