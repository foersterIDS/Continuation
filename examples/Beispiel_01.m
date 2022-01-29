clear;
close all;
clc;
%% move to folder of Beispiel_01.m
if(~isdeployed)
  cd(fileparts(which('Beispiel_01.m')));
end
addpath('..\src');
addpath('test_cases');
%% Testfunctions:
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
%% Solve:
[var_all,l_all,exitflag,bifs,s_all,last_jacobian,break_fun_out] = ...
    continuation(fun,v0,lams,lame,ds0,'bifurcation','mark','ds_max',ds_max,'plot','on','step_size_control','multiplicative');

%% Smoothen path
% [X_interp] = aux.smoothen_path([var_all; l_all], s_all);

%% plot smooth curves into figure
% plot_method = 'on'; % 'on' or 'three_dim'
% 
% l_interp = X_interp(end, :);
% var_interp = X_interp(1:end-1,:);
% 
% if ~strcmp(plot_method, 'three_dim')
%     num_pl = size(var_interp,1);
% else
%     num_pl = 1;
% end
% 
% colors = cell(num_pl,1);        
% for k = 1:num_pl
%     colors(k) = {plot.get_RGB(k,num_pl+5,5)};
% end
% 
% hold on;
% if strcmp(plot_method, 'three_dim')
%     pl = plot3(l_interp, var_interp(1,:), var_interp(2,:), 'LineStyle', ':', 'LineWidth', 2);
%     set(pl, {'Color'}, colors);
% else
%     pl = plot(l_interp, var_interp, 'LineStyle', ':', 'LineWidth', 2);
%     set(pl, {'Color'}, colors);
% end