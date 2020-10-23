clear all;
close all;
clc;
addpath('..\src');
addpath('test_cases');
%% Testfunktion #01
% 0=!v.^2+5-exp((1/l)*v)
% testfun01;
%% Testfunktion #02
% Duffing: mu \ddot q + zeta \dot q + kappa q + \gamma q^3 = P cos( Om * t )
testfun02;
%% Testfunktion #03
% Stochastic Duffing
% testfun03;
%% Testfunktion #04
% Point of intersection of circle and sin(radius)-scaled exponential function with radius as parameter
% testfun04;
%% Testfunktion #05
% function with bifurcation
% testfun05;
%% Testfunktion #06
% tests ellipsoid2
% testfun06;
%% Testfunktion #07
% Multidimensional (nd = 500)
% testfun07;
%% Solve
[vs,ls,exitflag,bifs] = continuation(fun,v0,lams,lame,ds0,'homotopy','on','solver','fsolve','bifurcation','mark','ds_max',ds_max,'plot','on');