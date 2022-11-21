clear;
close all;
clc;
addpath('../src/');
%% residual function:
fun = @(var,l) var-l;
% It is possible that fun returns R or [R,J] with R being the residual and
% J the Jacobian with respect to var or [var;l].
%% settings:
var0 = 0; % initial value of var
lStart = 0; % initial value of l
lEnd = 5; % end value of l
ds0 = 0.01; % initial step size
dsMin = 1*10^-5;
dsMax = 2*10^-1;

%% event based stepsize control with two simple events
% create StepSizeControlEvent object
eventObj = StepSizeControlEvent();

% Add first event
% set new dsMax to be 5e-2 when lambda is between 1 and 1.4
eventCondition = @(Path,dsCurrent) Path.lAll(end)+dsCurrent/3 >= 1 && Path.lAll(end)+dsCurrent/3 <= 1.4;
eventNeededParameters = {'Path','dsCurrent'};
eventDsMin = dsMin;
eventDsMax = 1e-3;
eventCounterMax = 1;
eventObj = eventObj.addEvent(eventCondition,eventNeededParameters,eventDsMin,eventDsMax,eventCounterMax);

% Add second event
% set new dsMax to be 1 when lambda is greater than 3
eventCondition = @(Path,dsCurrent) Path.lAll(end)+dsCurrent/3 > 3;
eventNeededParameters = {'Path','dsCurrent'};
eventDsMin = dsMin;
eventDsMax = 1e-2;
eventCounterMax = 1;
eventObj = eventObj.addEvent(eventCondition,eventNeededParameters,eventDsMin,eventDsMax,eventCounterMax);


%% run continuation:

% basic options & plot on (predictor tangential, dsMin/max) & activate
% event based step size control
[varAll,lAll,exitflag,~,sAll] = continuation(fun,var0,lStart,lEnd,ds0,'plot','on',...
                                                                      'predictor','tangential',...
                                                                      'dsMin',dsMin,'dsMax',dsMax,...
                                                                      'stepSizeEvent',true,...
                                                                      'eventUserInput',eventObj.getEvents);
%% Also plot used stepsize
figure('Units','normalized','Position',[0.2,0.2,0.6,0.6]);

ds = [0,sAll(2:end) - sAll(1:end-1)];
pl(1) = plot(lAll, ds,'bo',LineWidth=2); hold on;
xlabel('$\lambda$',Interpreter='latex');
ylabel('$\Delta s(\lambda)$',Interpreter='latex');
xl = xlim();
yl = ylim();
pl(2) = fill([1, 1.4, 1.4, 1],[yl(1),yl(1),yl(2),yl(2)],'r','EdgeColor','none','FaceAlpha',0.4); hold on;
pl(3) = fill([3, xl(2), xl(2), 3],[yl(1),yl(1),yl(2),yl(2)],'g','EdgeColor','none','FaceAlpha',0.4);
pl(4) = fill([xl(1), 1, 1, xl(1)],[yl(1),yl(1),yl(2),yl(2)],[0.5,0.5,0.5],'EdgeColor','none','FaceAlpha',0.4);
fill([1.4, 3, 3, 1.4],[yl(1),yl(1),yl(2),yl(2)],[0.5,0.5,0.5],'EdgeColor','none','FaceAlpha',0.4);
uistack(pl(1),'top');
legend(pl(1:4),'used stepsize','Event #1', 'Event #2', 'no event active');
fontsize(gcf, 13, 'points')