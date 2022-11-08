clear;
close all;
clc;
%% move to folder of Beispiel_03.m
if(~isdeployed)
  cd(fileparts(which('Beispiel_03.m')));
end
addpath('..\src');
addpath('test_cases');
%% Test functions:
testfun03; % Stochastic Duffing

%% Prepare stepsize event object:
% create StepSizeControlEvent-Object
event_obj = StepSizeControlEvent();

% Add first event
event_condition = @(Path) Path.l_all(end) >= 1 && Path.l_all(end) <= 1.1;
event_needed_parameters = {'Path'};
event_ds_min = 1e-10;
event_ds_max = 1e-3;
event_counter_max = 1;
event_obj = event_obj.addEvent(event_condition,event_needed_parameters,event_ds_min,event_ds_max,event_counter_max);

% Add second event
event_condition = @(Path) Path.l_all(end) >= 1.4;
event_needed_parameters = {'Path'};
event_ds_min = 1e-10;
event_ds_max = 1e-2;
event_counter_max = 1;
event_obj = event_obj.addEvent(event_condition,event_needed_parameters,event_ds_min,event_ds_max,event_counter_max);

%% Solve:
[var_all,l_all,exitflag,bifs,s_all,last_jacobian,break_fun_out] = ...
    continuation(fun,v0,lams,lame,ds0,'ds_max',ds_max,'plot','on','step_size_event',true,'event_user_input',event_obj.getEvents);