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
testfun12; % parabola intersecting lines
%% Solve:
[varAll,lAll,exitflag,bifs,sAll,jacobianOut,breakFunOut] = ...
    continuation(fun,v0,lams,lame,ds0,'dsMax',dsMax,'plot','off','corrector','sphere','bifurcation','trace');

%% Show results
fig = IDSfigure(1,'type','live'); clf;
hold on;
grid on; box on;
xlabel('$x$');
ylabel('$y$');
xlim([lams,lame]);

% Plot analytic solution for reference
lamPlot = linspace(lams,lame,1000);
parabolaPoints = sqrt(lamPlot);
idxReal = lamPlot>=0;
pl(1) = plot(lamPlot(idxReal),sqrt(lamPlot(idxReal)),'-','Color',getRGB('r'));
pl(2) = plot(lamPlot(idxReal),-sqrt(lamPlot(idxReal)),'-','Color',getRGB('r'));
pl(3) = plot(lamPlot,lamPlot/3,'-','Color',getRGB('r'));
pl(4) = plot(lamPlot,4-lamPlot,'-','Color',getRGB('r'));

% plot calulated path
plot(lAll,varAll(1,:),'-','Color',getRGB('b'));

% plot analytic bifs
r1 = roots([1,1,-4]).';
r2 = [0,3];

analyticBifs = [3, r1.^2, r2.^2;
                1, r1, r2];

for kk = 1:size(analyticBifs,2)
    plot(analyticBifs(1,kk),analyticBifs(2,kk),'x','Color',getRGB('r'),'MarkerSize',10);
end

% plot found bifs
foundBifs = bifs.bif;
for kk = 1:size(foundBifs,2)
    switch foundBifs(2,kk)
        case 0
            marker = 'x';
        case 1
            marker = 'o';
        otherwise
            marker = 's';
    end
    xBif = lAll(foundBifs(1,kk));
    yBif = varAll(1,foundBifs(1,kk));
    plot(xBif,yBif,'Color',getRGB('b'),'MarkerSize',10,'Marker',marker);
end

% legend
plL(1) = plot(NaN,NaN,'-','Color',getRGB('b'));
plL(2) = plot(NaN,NaN,'-','Color',getRGB('r'));

legend(plL,'found','missing');