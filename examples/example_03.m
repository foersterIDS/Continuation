clear;
close all;
clc;
addpath('../src/');
%% residual function:
fun = @(var,l) (l-2*var-1)*(l-(var-2).^2-5);
% It is possible that fun returns R or [R,J] with R being the residual and
% J the Jacobian with respect to var or [var;l].
%% settings:
var0 = 1; % initial value of var
lStart = 0; % initial value of l
lEnd = 10; % end value of l
ds0 = 0.01; % initial step size
dsMax = 0.1; % max. step size
%% run continuation:

% show plot and use dsMax
[varAll,lAll,exitflag] = continuation(fun,var0,lStart,lEnd,ds0,'plot','on',...
                                                                   'dsMax',dsMax);
input('Press any key to continue...');

% show plot, use dsMax and turn bifurcation detection off
[varAll,lAll,exitflag] = continuation(fun,var0,lStart,lEnd,ds0,'plot','on',...
                                                                   'dsMax',dsMax,...
                                                                   'bifurcation','off');
input('Press any key to continue...');

% show plot, use dsMax and set bifurcation detection to mark
[varAll,lAll,exitflag] = continuation(fun,var0,lStart,lEnd,ds0,'plot','on',...
                                                                   'dsMax',dsMax,...
                                                                   'bifurcation','mark');
input('Press any key to continue...');

% show plot, use dsMax and set bifurcation detection to determine
[varAll,lAll,exitflag] = continuation(fun,var0,lStart,lEnd,ds0,'plot','on',...
                                                                   'dsMax',dsMax,...
                                                                   'bifurcation','determine');
input('Press any key to continue...');

% show plot, use dsMax and set bifurcation detection to trace
[varAll,lAll,exitflag] = continuation(fun,var0,lStart,lEnd,ds0,'plot','on',...
                                                                   'dsMax',dsMax,...
                                                                   'bifurcation','trace');
input('Press any key to continue...');

%% dpa residual function:
funDpa = @(var,l,g) (l-2*var-g)*(l-(var-2*g).^2-5);
% It is possible that fun returns R or [R,J] with R being the residual and
% J the Jacobian with respect to var or [var;l].
%% dpa settings:
lEnd = 25; % end value of l
g0 = 1; % initial value of g
gTarget = 3; % end value of g
ds0 = 0.001; % initial step size
dsMax = 0.05; % max. step size
%% run continuation:

% show dpa plot, use dsMax and set bifurcation detection to parameterTrace
[varAllDpa,lAllDpa,exitflagDpa] = continuation(funDpa,var0,lStart,lEnd,ds0,'plot','dpa',...
                                                                           'dsMax',dsMax,...
                                                                           'bifurcation','parameterTrace',...
                                                                           'g0',g0,...
                                                                           'gTarget',gTarget);

% show plot, use dsMax and set bifurcation detection to trace, fix at g0
[varAllTraceG0,lAllTraceG0,exitflagTraceG0] = continuation(@(v,l) funDpa(v,l,g0),var0,lStart,lEnd,ds0,'plot','off',...
                                                                                                      'dsMax',dsMax,...
                                                                                                      'bifurcation','trace');

% show plot, use dsMax and set bifurcation detection to trace, fix at gTarget
[varAllTraceGTarget,lAllTraceGTarget,exitflagTraceGTarget] = continuation(@(v,l) funDpa(v,l,gTarget),var0,lStart,lEnd,ds0,'plot','off',...
                                                                                                                          'dsMax',dsMax,...
                                                                                                                          'bifurcation','trace');

%% plot:
figure(8);
clf;
plot3(lAllDpa(1,:),lAllDpa(2,:),varAllDpa,'k-','LineWidth',2);
hold on;
plot3(lAllTraceG0,g0*ones(size(lAllTraceG0)),varAllTraceG0,'b-','LineWidth',2);
plot3(lAllTraceGTarget,gTarget*ones(size(lAllTraceGTarget)),varAllTraceGTarget,'r-','LineWidth',2);
hold off;
grid on;
xlim([lStart,lEnd]);
ylim([g0,gTarget]);
drawnow;