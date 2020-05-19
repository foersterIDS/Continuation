clear all;
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
% Stochastic (Duffing or ECS)
% testfun03;
%% Solve
[vs,ls,exitflag] = continuation(fun,v0,lams,lame,ds0,'homotopy','on','solver','newton');%,'arclength','ellipsoid');
%% Plot
figure(1);
clf;
plot(ls,vs,'-','LineWidth',2);
hold off;
grid on;
xlim([lams,lame]);
xlabel('$\lambda$','interpreter','latex');
ylabel('$v_{i}$','interpreter','latex');
drawnow;