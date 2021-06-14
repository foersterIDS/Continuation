clear all;
close all;
clc;
addpath('..\src');
addpath('test_cases');
%% Testfunktions:
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
%% Solve:
[var_all,l_all,exitflag,bifs,s_all,last_jacobian,break_fun_out] = ...
    continuation(fun,v0,lams,lame,ds0,'bifurcation','trace','ds_max',ds_max,'plot','detail', 'plot_pause', 10);