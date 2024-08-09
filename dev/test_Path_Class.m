clear;
close all;
clc;
addpath('../src/')
%% run class test
Opt = struct('jacobianOut',struct('full',true,'basic',false)); % full/basic
OptIsSet = struct('pathInfoFunction',true,...
                  'bifAdditionalTestfunction',true);
Solver = struct('output',struct('iterations',round(1+10*rand),'rateOfContraction',rand));
StepsizeOptions = struct('iterations',true,...
                         'predictor',true,...
                         'rateOfContraction',true,...
                         'speedOfContinuation',true);

cPath = continuation.Path(3,1,Opt,OptIsSet,StepsizeOptions);

cPath.addPointAtEnd(randn(3,1),1,randn(3,4),Solver,'predictor',randn(3,1),'speedOfContinuation',rand,'pathInfoValue',rand,'bifTestValue',rand);
cPath.addPointAtEnd(randn(3,1),2,randn(3,4),Solver,'predictor',randn(3,1),'speedOfContinuation',rand,'pathInfoValue',rand,'bifTestValue',rand);
cPath.addPointAtEnd(randn(3,1),3,randn(3,4),Solver,'predictor',randn(3,1),'speedOfContinuation',rand,'pathInfoValue',rand,'bifTestValue',rand);
cPath.toggleStepback();
cPath.addPointAtEnd(randn(3,1),4,randn(3,4),Solver,'predictor',randn(3,1),'speedOfContinuation',rand,'pathInfoValue',rand,'bifTestValue',rand);
cPath.toggleStepback();
cPath.addPointAtEnd(randn(3,1),5,randn(3,4),Solver,'predictor',randn(3,1),'speedOfContinuation',rand,'pathInfoValue',rand,'bifTestValue',rand);
cPath.remove(5);
cPath.addPointAtEnd(randn(3,1),6,randn(3,4),Solver,'predictor',randn(3,1),'speedOfContinuation',rand,'pathInfoValue',rand,'bifTestValue',rand);
cPath.addPointAtEnd(randn(3,1),7,randn(3,4),Solver,'predictor',randn(3,1),'speedOfContinuation',rand,'pathInfoValue',rand,'bifTestValue',rand);
cPath.remove(5:6);
cPath.addPointAtEnd(randn(3,1),8,randn(3,4),Solver,'predictor',randn(3,1),'speedOfContinuation',rand,'pathInfoValue',rand,'bifTestValue',rand);
cPath.addPointAtEnd(randn(3,1),9,randn(3,4),Solver,'predictor',randn(3,1),'speedOfContinuation',rand,'pathInfoValue',rand,'bifTestValue',rand);
cPath.addPointAtEnd(randn(3,1),10,randn(3,4),Solver,'predictor',randn(3,1),'speedOfContinuation',rand,'pathInfoValue',rand,'bifTestValue',rand);
cPath.addPointAtEnd(randn(3,1),11,randn(3,4),Solver,'predictor',randn(3,1),'speedOfContinuation',rand,'pathInfoValue',rand,'bifTestValue',rand);
cPath.remove(1:2);
cPath.addPoint(randn(3,1),1,randn(3,4),Solver,1,'predictor',randn(3,1),'speedOfContinuation',rand,'pathInfoValue',rand,'bifTestValue',rand);
cPath.addPoint(randn(3,1),2,randn(3,4),Solver,2,'predictor',randn(3,1),'speedOfContinuation',rand,'pathInfoValue',rand,'bifTestValue',rand);
cPath.remove(3);
cPath.addPoint(randn(3,1),4,randn(3,4),Solver,4,'predictor',randn(3,1),'speedOfContinuation',rand,'pathInfoValue',rand,'bifTestValue',rand);
cPath.addPoint(randn(3,1),5,randn(3,4),Solver,5,'predictor',randn(3,1),'speedOfContinuation',rand,'pathInfoValue',rand,'bifTestValue',rand);
cPath.addPoint(randn(3,1),6,randn(3,4),Solver,6,'predictor',randn(3,1),'speedOfContinuation',rand,'pathInfoValue',rand,'bifTestValue',rand);
cPath.addPoint(randn(3,1),7,randn(3,4),Solver,7,'predictor',randn(3,1),'speedOfContinuation',rand,'pathInfoValue',rand,'bifTestValue',rand);
cPath.remove(11);

% 1:10
cPath.lAll

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