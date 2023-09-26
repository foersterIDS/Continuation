clear;
close all;
clc;
%% move to folder of Beispiel_01.m
if(~isdeployed)
  cd(fileparts(which('Beispiel_01.m')));
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
% testfun11; % basic test solverForce1it
% testfun12; % parabola intersecting lines
% testfun13; % y = sin(1/x)
% testfun14; % Infinity, Test for (...,'plot','threeDim',...)
%% Solve:
[varAll,lAll,exitflag,bifs,sAll,jacobianOut,breakFunOut] = ...
    continuation(fun,v0,lams,lame,ds0,'dsMax',dsMax,'plot','on','corrector','sphere');

% %% Smoothen path
% XInterp = aux.smoothenPath([varAll; lAll], sAll);
% 
% % plot smooth curves into figure
% plotMethod = 'on'; % 'on' or 'threeDim'
% 
% axPrev = gca;
% 
% fig = figure('units', 'normalized', 'position', [0.2,0.3,0.6,0.5]);
% 
% copyobj(axPrev,fig);
% 
% lInterp = XInterp(end, :);
% varInterp = XInterp(1:end-1,:);
% 
% if ~strcmp(plotMethod, 'threeDim')
%     numPl = size(varInterp,1);
% else
%     numPl = 1;
% end
% 
% colors = cell(numPl,1);        
% for k = 1:numPl
%     colors(k) = {plot.getRGB(k,numPl,1)};
% end
% 
% hold on;
% cla;
% if strcmp(plotMethod, 'threeDim')
%     pl = plot3(lInterp, varInterp(1,:), varInterp(2,:), 'LineStyle', '-', 'LineWidth', 2);
%     set(pl, {'Color'}, colors);
% else
%     pl = plot(lInterp, varInterp, 'LineStyle', '-', 'LineWidth', 2);
%     set(pl, {'Color'}, colors);
% end