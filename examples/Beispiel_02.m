clear;
close all;
clc;
%% move to folder of Beispiel_02.m
if(~isdeployed)
  cd(fileparts(which('Beispiel_02.m')));
end
addpath('../src');
addpath('testCases');
%% Test functions:
% testfun01; % 0=!v.^2+5-exp((1/l)*v)
testfun02; % Duffing: mu \ddot q + zeta \dot q + kappa q + \gamma q^3 = P cos( Om * t )
% testfun03; % Stochastic Duffing
% testfun04; % Point of intersection of circle and sin(radius)-scaled exponential function with radius as parameter
% testfun05; % function with bifurcation
% testfun06; % tests ellipsoid2
% testfun07; % Multidimensional (nd = 500)
% testfun08; % Pitchfork-Bifurkation
% testfun09; % circle
% testfun10; % NB-excitation gap-oszillator
% testfun11; % basic test solver_force1it
% testfun12; % parabola intersecting lines
% testfun13; % y = sin(1/x)
% testfun14; % Infinity ('plot','three_dim')
%% Prepare for homotopy:
ll = abs(lame-lams)/2+min(lams,lame); % aimed value of lambda
funh = @(var) fun(var,ll);
%% Solve:
[v_solution,exitflag] = homotopy(funh,v0,'plot','on','homotopy','fix');