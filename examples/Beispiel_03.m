clear;
close all;
clc;
%% move to folder of Beispiel_03.m
if(~isdeployed)
  cd(fileparts(which('Beispiel_03.m')));
end
addpath('..\src');
addpath('testCases');
%% Test functions:
testfun03; % Stochastic Duffing

%% Prepare stepsize event object:
% create StepSizeControlEvent-Object
eventObj = StepSizeControlEvent();

% Add first event
eventCondition = @(Path) Path.lAll(end) >= 1 && Path.lAll(end) <= 1.1;
eventNeededParameters = {'Path'};
eventDsMin = 1e-10;
eventDsMax = 1e-3;
eventCounterMax = 1;
eventObj = eventObj.addEvent(eventCondition,eventNeededParameters,eventDsMin,eventDsMax,eventCounterMax);

% Add second event
eventCondition = @(Path) Path.lAll(end) >= 1.4;
eventNeededParameters = {'Path'};
eventDsMin = 1e-10;
eventDsMax = 1e-2;
eventCounterMax = 1;
eventObj = eventObj.addEvent(eventCondition,eventNeededParameters,eventDsMin,eventDsMax,eventCounterMax);

%% Solve:
[varAll,lAll,exitflag,bifs,sAll,lastJacobian,breakFunOut] = ...
    continuation(fun,v0,lams,lame,ds0,'dsMax',dsMax,'plot','on','stepSizeEvent',true,'eventUserInput',eventObj.getEvents);