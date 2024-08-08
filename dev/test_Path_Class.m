clear;
close all;
clc;
addpath('../src/')
%% run class test
Opt = struct('jacobianOut','full'); % full/basic
OptIsSet = struct('pathInfoFunction',false,...
                  'bifAdditionalTestfunction',false);
StepsizeOptions = struct('iterations',true,...
                         'predictor',true,...
                         'rateOfContraction',true,...
                         'speedOfContinuation',true);

p = continuation.Path(3,1,Opt,OptIsSet,StepsizeOptions);

p.addPointAtEnd(randn(3,1),1,randn(3,4),'iterations',round(1+10*rand),'predictor',randn(3,1),'rateOfContraction',rand,'speedOfContinuation',rand);
p.addPointAtEnd(randn(3,1),2,randn(3,4),'iterations',round(1+10*rand),'predictor',randn(3,1),'rateOfContraction',rand,'speedOfContinuation',rand);
p.addPointAtEnd(randn(3,1),3,randn(3,4),'iterations',round(1+10*rand),'predictor',randn(3,1),'rateOfContraction',rand,'speedOfContinuation',rand);
p.toggleStepback();
p.addPointAtEnd(randn(3,1),4,randn(3,4),'iterations',round(1+10*rand),'predictor',randn(3,1),'rateOfContraction',rand,'speedOfContinuation',rand);
p.toggleStepback();
p.addPointAtEnd(randn(3,1),5,randn(3,4),'iterations',round(1+10*rand),'predictor',randn(3,1),'rateOfContraction',rand,'speedOfContinuation',rand);
p.remove(5);
p.addPointAtEnd(randn(3,1),6,randn(3,4),'iterations',round(1+10*rand),'predictor',randn(3,1),'rateOfContraction',rand,'speedOfContinuation',rand);
p.addPointAtEnd(randn(3,1),7,randn(3,4),'iterations',round(1+10*rand),'predictor',randn(3,1),'rateOfContraction',rand,'speedOfContinuation',rand);
p.remove(5:6);
p.addPointAtEnd(randn(3,1),8,randn(3,4),'iterations',round(1+10*rand),'predictor',randn(3,1),'rateOfContraction',rand,'speedOfContinuation',rand);
p.addPointAtEnd(randn(3,1),9,randn(3,4),'iterations',round(1+10*rand),'predictor',randn(3,1),'rateOfContraction',rand,'speedOfContinuation',rand);
p.addPointAtEnd(randn(3,1),10,randn(3,4),'iterations',round(1+10*rand),'predictor',randn(3,1),'rateOfContraction',rand,'speedOfContinuation',rand);
p.addPointAtEnd(randn(3,1),11,randn(3,4),'iterations',round(1+10*rand),'predictor',randn(3,1),'rateOfContraction',rand,'speedOfContinuation',rand);
p.remove(1:2);

% 4,3,8,9,10,11

p.xAll

p.xPredictor = p.xAll(:,end);

%% run continuation test
% add paths
addpath('../src');
addpath('../examples/testCases/');
addpath('../examples/testCases/testfun02_aux/');
% Test functions:
testfun02; % Duffing: mu \ddot q + zeta \dot q + kappa q + \gamma q^3 = P cos( Om * t )
% Solve:
[varAll,lAll,exitflag,bifs,sAll,jacobianOut,breakFunOut] = ...
    continuation(fun,v0,lams,lame,ds0,'dsMax',dsMax,'plot','on','corrector','sphere','pauseOnError',true);