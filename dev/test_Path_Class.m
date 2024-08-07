clear;
clc;
%% run class test
p = continuation.Path(3,1,"saveAllJacobian",true);

p.addPointAtEnd(randn(3,1),randn,randn(3,4));
p.addPointAtEnd(randn(3,1),randn,randn(3,4));
p.addPointAtEnd(randn(3,1),randn,randn(3,4));
p.addPointAtEnd(randn(3,1),randn,randn(3,4));
p.addPointAtEnd(randn(3,1),randn,randn(3,4));
p.addPointAtEnd(randn(3,1),randn,randn(3,4));
p.addPointAtEnd(randn(3,1),randn,randn(3,4));
p.addPointAtEnd(randn(3,1),randn,randn(3,4));

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