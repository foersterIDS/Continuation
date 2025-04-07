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
testfun02;
%% Solve:
%% figure vorab erstellen und handle als globale Variable speichern
fHandle = figure(1);
fHandle.Color = 'w';
fHandle.Units = 'centimeters';
fHandle.Position(3:4) = [24.8, 11.56];
axis([lams,lame,-6,9]); grid on; box on;
xlabel('$x$','Interpreter','latex'); ylabel('$y$','Interpreter','latex');

%% Optionen für Animation festlegen
% Anzahl der Pfadverfolgungspunkte, die zwischendurch nicht exportiert
% werden soll (=0, dann werden alle Aktualisierungen dargestellt)
%
skipFrames = 5; 
%
% Abweichende DPI, (=NaN, dann wird die Auflösung des Bildschirms
% verwendet)
%
customDPI = 200;
%
% Name des Videos
%
fName = 'Test';
%
% FPS
%
fps = 30;

%% Pfadverfolgung aufrufen
[varAll,lAll,exitflag,bifs,sAll,jacobianOut,breakFunOut] = ...
    continuation(fun,v0,lams,lame,ds0,'dsMax',dsMax,'plot','on','corrector','sphere',...
    'plotOptions',plot.PlotOptions("figure",fHandle,'createAnimation',true,'skipFrames',skipFrames,'customDPI',customDPI,'animationFilename',fName,'fps',fps));