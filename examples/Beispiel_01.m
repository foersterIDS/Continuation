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
%% Testfunktion #04
% Point of intersection of circle and sin(radius)-scaled exponential function with radius as parameter
% testfun04;
%% Testfunktion #05
% function with bifurcation
% testfun05;
%% Solve
[vs,ls,exitflag,bifs] = continuation(fun,v0,lams,lame,ds0,'homotopy','on','solver','newton');
%% Plot
figure(1);
clf;
plot(ls,vs,'-','LineWidth',2);
if ~isempty(bifs)
    hold on;
    plot(ls(bifs(1,bifs(2,:)==0)),vs(:,bifs(1,bifs(2,:)==0)),'bo','LineWidth',2);
    plot(ls(bifs(1,bifs(2,:)==1)),vs(:,bifs(1,bifs(2,:)==1)),'rx','LineWidth',2);
    hold off;
end
grid on;
xlim([lams,lame]);
xlabel('$\lambda$','interpreter','latex');
ylabel('$v_{i}$','interpreter','latex');
drawnow;