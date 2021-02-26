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
% testfun11; % basic test
%% Solve:
[vs,ls,exitflag,bifs] = continuation(fun,v0,lams,lame,ds0,'bifurcation','off','ds_max',ds_max,'plot','on');