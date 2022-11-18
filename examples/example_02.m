clear;
close all;
clc;
addpath('../src/');
%% residual function:
fun = @(var,l) var.^2+l.^2-1;
% It is possible that fun returns R or [R,J] with R being the residual and
% J the Jacobian with respect to var or [var;l].
%% settings:
var0 = 1/sqrt(2); % initial value of var
lStart = -2; % left boundary of l
lEnd = +2; % right boundary of l
ds0 = 0.001; % initial step size
l0 = 1/sqrt(2); % initial value of l
%% run continuation:

% start at l0 and show plot
[varAll,lAll,exitflag] = continuation(fun,var0,lStart,lEnd,ds0,'l0',l0,...
                                                                   'plot','on');
input('Press any key to continue...');

% start at l0 and show plot, turn closedCurveDetection off and
% set nStepMax to 500
[varAll,lAll,exitflag] = continuation(fun,var0,lStart,lEnd,ds0,'l0',l0,...
                                                                   'plot','on',...
                                                                   'closedCurveDetection','off',...
                                                                   'nStepMax',500);