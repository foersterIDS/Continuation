%% continuation logo
clear; close all; clc;

h = figure();
axes('Units', 'normalized', 'Position', [0 0 1 1]);
t = -2.6:0.5:2.4;
x = 0.5*t.^3-2*t;

plot(t,x,'b','LineWidth',4); hold on;
plot(t,x,'ro','MarkerFaceColor','r', 'MarkerSize',14);
axis([-3,3,-8,8]);
set(gca,'visible','off');

h.Color = [1,1,1];