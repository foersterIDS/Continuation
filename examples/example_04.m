clear;
close all;
clc;
addpath('../src/');
%% residual function:
fun = @(var,l) var-l*sin(10/(l-2.5));
% It is possible that fun returns R or [R,J] with R being the residual and
% J the Jacobian with respect to var or [var;l].
%% settings:
var0 = 0; % initial value of var
lStart = 0; % initial value of l
lEnd = 2; % end value of l
ds0 = 0.01; % initial step size
dsMin = 1*10^-5;
dsMax = 2*10^-1;
ssc = {'angleChange';
       'angleCustom';
       'contraction';
       'error';
       'errorAlt';
       'fayezioghani';
       'fix';
       'iterationsExponential';
       'iterationsPolynomial';
       'multiplicative'; % inital value
       'multiplicativeAlt';
       'pidCustom';
       'pidValli';
       'szyszkowski';
       'yoon'};
nSsc = numel(ssc);
alphaReverse = 2.5;
%% run continuation:

% basic options & plot on (alphaReverse pi, predictor tangential, dsMin/max)
[varAll0,lAll0,exitflag0] = continuation(fun,var0,lStart,lEnd,ds0,'plot','on',...
                                                                         'alphaReverse',alphaReverse,...
                                                                         'predictor','tangential',...
                                                                         'dsMin',dsMin,'dsMax',dsMax);
drawnow;

% loop over all step size control methods
clc;
figure(2);
clf;
ax = axis;
xlim([lStart,lEnd]);
ylim([-2,2]);
xlabel('$\lambda$','interpreter','latex');
ylabel('$v$','interpreter','latex');
grid on;
box on;
drawnow;
for ii=1:nSsc
    fprintf('run with scc: %s',ssc{ii});
    tic;
    [varAll{ii},lAll{ii},exitflag(ii)] = continuation(fun,var0,lStart,lEnd,ds0,'plot','off',...
                                                                                   'display','off',...
                                                                                   'alphaReverse',alphaReverse,...
                                                                                   'predictor','tangential',...
                                                                                   'dsMin',dsMin,'dsMax',dsMax,...
                                                                                   'stepSizeControl',ssc{ii});
    tSsc(ii) = toc;
    fprintf(' (duration: %.3f s, target reached: %d)\n',tSsc(ii),lAll{ii}(end)>=lEnd);
    hold on;
    plot(lAll{ii},varAll{ii},'-','LineWidth',2,'Color',plot.getRGB(ii,nSsc,1));
    hold off;
    drawnow;
end

% result:
[tSscMin,iiMin] = min(tSsc);
fprintf('- - - - - - - - - - - - - - - - - - - - - - - - - - - -\nFastest method: %s (%.3f s)\n',ssc{iiMin},tSscMin);