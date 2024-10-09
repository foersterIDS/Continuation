clear;
close all;
clc;
tic
addpath('..\src');
addpath('testCases');
addpath('crashtestAux');
probinfo = [];
probcounter = 0;
%% Test function #04
% Point of intersection of circle and sin(radius)-scaled exponential function with radius as parameter
% testfun04;
%% Test function #08
% Pitchfork-Bifurkation
testfun08;
%% Crashtest:
for i=1:2
    % i=1: with jacobian, i=2: without jacobian
    %% standard-config.:
    fprintf('\n### %d: standard-config ###\n',i);
    [vs,ls,exitflag] = continuation(funJacoTest{i},v0,lams,lame,ds0,'dsMax',dsMax);
    probinfo = crashtestCheckOutput(vs,ls,exitflag,lams,lame,sprintf('\n### %d: standard-config ###\n',i),probinfo,probcounter);
    %% check residual on:
    fprintf('\n### %d: check residual: on ###\n',i);
    [vs,ls,exitflag] = continuation(funJacoTest{i},v0,lams,lame,ds0,'dsMax',dsMax,'checkResidual','on');
    [probinfo,probcounter] = crashtestCheckOutput(vs,ls,exitflag,lams,lame,sprintf('\n### %d: check residual: on ###\n',i),probinfo,probcounter);
    %% display off:
    fprintf('\n### %d: correct predictor: off ###\n',i);
    [vs,ls,exitflag] = continuation(funJacoTest{i},v0,lams,lame,ds0,'dsMax',dsMax,'correctPredictor','off');
    [probinfo,probcounter] = crashtestCheckOutput(vs,ls,exitflag,lams,lame,sprintf('\n### %d: correct predictor: off ###\n',i),probinfo,probcounter);
    %% display off:
    fprintf('\n### %d: display: off ###\n',i);
    [vs,ls,exitflag] = continuation(funJacoTest{i},v0,lams,lame,ds0,'dsMax',dsMax,'display','off');
    [probinfo,probcounter] = crashtestCheckOutput(vs,ls,exitflag,lams,lame,sprintf('\n### %d: display: off ###\n',i),probinfo,probcounter);
    %% homotopy:
    fprintf('\n### %d: homotopy: off ###\n',i);
    [vs,ls,exitflag] = continuation(funJacoTest{i},v0,lams,lame,ds0,'dsMax',dsMax,'display','off','homotopy','off');
    [probinfo,probcounter] = crashtestCheckOutput(vs,ls,exitflag,lams,lame,sprintf('\n### %d: homotopy: off ###\n',i),probinfo,probcounter);
    fprintf('\n### %d: homotopy: on ###\n',i);
    [vs,ls,exitflag] = continuation(funJacoTest{i},v0,lams,lame,ds0,'dsMax',dsMax,'display','off','homotopy','on');
    [probinfo,probcounter] = crashtestCheckOutput(vs,ls,exitflag,lams,lame,sprintf('\n### %d: homotopy: on ###\n',i),probinfo,probcounter);
    fprintf('\n### %d: homotopy: fix ###\n',i);
    [vs,ls,exitflag] = continuation(funJacoTest{i},v0,lams,lame,ds0,'dsMax',dsMax,'display','off','homotopy','fix');
    [probinfo,probcounter] = crashtestCheckOutput(vs,ls,exitflag,lams,lame,sprintf('\n### %d: homotopy: fix ###\n',i),probinfo,probcounter);
    fprintf('\n### %d: homotopy: fixnt ###\n',i);
    [vs,ls,exitflag] = continuation(funJacoTest{i},v0,lams,lame,ds0,'dsMax',dsMax,'display','off','homotopy','fixnt');
    [probinfo,probcounter] = crashtestCheckOutput(vs,ls,exitflag,lams,lame,sprintf('\n### %d: homotopy: fixnt ###\n',i),probinfo,probcounter);
    fprintf('\n### %d: homotopy: newton ###\n',i);
    [vs,ls,exitflag] = continuation(funJacoTest{i},v0,lams,lame,ds0,'dsMax',dsMax,'display','off','homotopy','newton');
    [probinfo,probcounter] = crashtestCheckOutput(vs,ls,exitflag,lams,lame,sprintf('\n### %d: homotopy: newton ###\n',i),probinfo,probcounter);
    fprintf('\n### %d: homotopy: f2 ###\n',i);
    [vs,ls,exitflag] = continuation(funJacoTest{i},v0,lams,lame,ds0,'dsMax',dsMax,'display','off','homotopy','f2');
    [probinfo,probcounter] = crashtestCheckOutput(vs,ls,exitflag,lams,lame,sprintf('\n### %d: homotopy: newton ###\n',i),probinfo,probcounter);
    %% solver:
    fprintf('\n### %d: solver: fsolve ###\n',i);
    [vs,ls,exitflag] = continuation(funJacoTest{i},v0,lams,lame,ds0,'dsMax',dsMax,'display','off','solver','fsolve');
    [probinfo,probcounter] = crashtestCheckOutput(vs,ls,exitflag,lams,lame,sprintf('\n### %d: solver: fsolve ###\n',i),probinfo,probcounter);
    fprintf('\n### %d: solver: lsqnonlin ###\n',i);
    [vs,ls,exitflag] = continuation(funJacoTest{i},v0,lams,lame,ds0,'dsMax',dsMax,'display','on','solver','lsqnonlin');
    [probinfo,probcounter] = crashtestCheckOutput(vs,ls,exitflag,lams,lame,sprintf('\n### %d: solver: lsqnonlin ###\n',i),probinfo,probcounter);
    fprintf('\n### %d: solver: newton ###\n',i);
    [vs,ls,exitflag] = continuation(funJacoTest{i},v0,lams,lame,ds0,'dsMax',dsMax,'display','off','solver','newton');
    [probinfo,probcounter] = crashtestCheckOutput(vs,ls,exitflag,lams,lame,sprintf('\n### %d: solver: newton ###\n',i),probinfo,probcounter);
    %% breakFunction:
    fprintf('\n### %d: breakFunction ###\n',i);
    [vs,ls,exitflag] = continuation(funJacoTest{i},v0,lams,lame,ds0,'dsMax',dsMax,'display','off','breakFunction',@(f,J,v,l,breakFunOut) aux.initialBreakFunction(f,J,v,l,breakFunOut));
    [probinfo,probcounter] = crashtestCheckOutput(vs,ls,exitflag,lams,lame,sprintf('\n### %d: breakFunction ###\n',i),probinfo,probcounter);
    %% corrector:
    fprintf('\n### %d: corrector: sphere ###\n',i);
    [vs,ls,exitflag] = continuation(funJacoTest{i},v0,lams,lame,ds0,'dsMax',dsMax,'display','off','corrector','sphere');
    [probinfo,probcounter] = crashtestCheckOutput(vs,ls,exitflag,lams,lame,sprintf('\n### %d: corrector: sphere ###\n',i),probinfo,probcounter);
    fprintf('\n### %d: corrector: orthogonal ###\n',i);
    [vs,ls,exitflag] = continuation(funJacoTest{i},v0,lams,lame,ds0,'dsMax',dsMax,'display','off','corrector','orthogonal');
    [probinfo,probcounter] = crashtestCheckOutput(vs,ls,exitflag,lams,lame,sprintf('\n### %d: corrector: orthogonal ###\n',i),probinfo,probcounter);
    fprintf('\n### %d: corrector: orthogonal2 ###\n',i);
    [vs,ls,exitflag] = continuation(funJacoTest{i},v0,lams,lame,ds0,'dsMax',dsMax,'display','off','corrector','orthogonal2');
    [probinfo,probcounter] = crashtestCheckOutput(vs,ls,exitflag,lams,lame,sprintf('\n### %d: corrector: orthogonal2 ###\n',i),probinfo,probcounter);
    fprintf('\n### %d: corrector: ellipsoid ###\n',i);
    [vs,ls,exitflag] = continuation(funJacoTest{i},v0,lams,lame,ds0,'dsMax',dsMax,'dsMin',10^-4,'display','on','corrector','ellipsoid','predictor','tangential');
    [probinfo,probcounter] = crashtestCheckOutput(vs,ls,exitflag,lams,lame,sprintf('\n### %d: corrector: ellipsoid ###\n',i),probinfo,probcounter);
    fprintf('\n### %d: corrector: ellipsoid2 ###\n',i);
    [vs,ls,exitflag] = continuation(funJacoTest{i},v0,lams,lame,ds0,'dsMax',dsMax,'nStepMax',1000,'display','off','corrector','ellipsoid2');
    [probinfo,probcounter] = crashtestCheckOutput(vs,ls,exitflag,lams,lame,sprintf('\n### %d: corrector: ellipsoid2 ###\n',i),probinfo,probcounter);
    fprintf('\n### %d: corrector: unique ###\n',i);
    [vs,ls,exitflag] = continuation(funJacoTest{i},v0,lams,lame,ds0,'dsMax',dsMax,'display','off','corrector','unique');
    [probinfo,probcounter] = crashtestCheckOutput(vs,ls,exitflag,lams,lame,sprintf('\n### %d: corrector: unique ###\n',i),probinfo,probcounter);
    fprintf('\n### %d: corrector: paraboloid ###\n',i);
    [vs,ls,exitflag] = continuation(funJacoTest{i},v0,lams,lame,ds0,'dsMax',dsMax,'display','off','corrector','paraboloid');
    [probinfo,probcounter] = crashtestCheckOutput(vs,ls,exitflag,lams,lame,sprintf('\n### %d: corrector: paraboloid ###\n',i),probinfo,probcounter);
    %% deflation:
    fprintf('\n### %d: deflation: on ###\n',i);
    [vs,ls,exitflag] = continuation(funJacoTest{i},v0,lams,lame,ds0,'dsMax',dsMax,'display','off','deflation','on');
    [probinfo,probcounter] = crashtestCheckOutput(vs,ls,exitflag,lams,lame,sprintf('\n### %d: deflation: on ###\n',i),probinfo,probcounter);
    fprintf('\n### %d: deflation: off ###\n',i);
    [vs,ls,exitflag] = continuation(funJacoTest{i},v0,lams,lame,ds0,'dsMax',dsMax,'display','off','deflation','off');
    [probinfo,probcounter] = crashtestCheckOutput(vs,ls,exitflag,lams,lame,sprintf('\n### %d: deflation: off ###\n',i),probinfo,probcounter);
    %% bifurcation:
    fprintf('\n### %d: bifurcation: mark ###\n',i);
    [vs,ls,exitflag] = continuation(funJacoTest{i},v0,lams,lame,ds0,'dsMax',dsMax,'display','off','bifurcation','mark');
    [probinfo,probcounter] = crashtestCheckOutput(vs,ls,exitflag,lams,lame,sprintf('\n### %d: bifurcation: mark ###\n',i),probinfo,probcounter);
    fprintf('\n### %d: bifurcation: determine ###\n',i);
    [vs,ls,exitflag] = continuation(funJacoTest{i},v0,lams,lame,ds0,'dsMax',dsMax,'display','off','bifurcation','determine');
    [probinfo,probcounter] = crashtestCheckOutput(vs,ls,exitflag,lams,lame,sprintf('\n### %d: bifurcation: determine ###\n',i),probinfo,probcounter);
    fprintf('\n### %d: bifurcation: trace ###\n',i);
    [vs,ls,exitflag] = continuation(funJacoTest{i},v0,lams,lame,ds0,'dsMax',dsMax,'display','off','bifurcation','trace','nStepMax',200);
    [probinfo,probcounter] = crashtestCheckOutput(vs,ls,exitflag,lams,lame,sprintf('\n### %d: bifurcation: trace ###\n',i),probinfo,probcounter);
    fprintf('\n### %d: stopOnBifurcation: on ###\n',i);
    [vs,ls,exitflag] = continuation(funJacoTest{i},v0,lams,lame,ds0,'dsMax',dsMax,'display','off','stopOnBifurcation','on');
    [probinfo,probcounter] = crashtestCheckOutput(vs,ls,exitflag,lams,lame,sprintf('\n### %d: stopOnBifurcation: on ###\n',i),probinfo,probcounter);
    %% nIterOpt:
    fprintf('\n### %d: nIterOpt: 5 ###\n',i);
    [vs,ls,exitflag] = continuation(funJacoTest{i},v0,lams,lame,ds0,'dsMax',dsMax,'display','off','nIterOpt',5);
    [probinfo,probcounter] = crashtestCheckOutput(vs,ls,exitflag,lams,lame,sprintf('\n### %d: nIterOpt: 5 ###\n',i),probinfo,probcounter);
    %% dsMin:
    fprintf('\n### %d: dsMin: ds0/2 ###\n',i);
    [vs,ls,exitflag] = continuation(funJacoTest{i},v0,lams,lame,ds0,'dsMax',dsMax,'display','off','dsMin',ds0/2);
    [probinfo,probcounter] = crashtestCheckOutput(vs,ls,exitflag,lams,lame,sprintf('\n### %d: dsMin: ds0/2 ###\n',i),probinfo,probcounter);
    %% dsTol:
    fprintf('\n### %d: dsTol: [0.5,1.5] ###\n',i);
    [vs,ls,exitflag] = continuation(funJacoTest{i},v0,lams,lame,ds0,'dsMax',dsMax,'display','off','dsTol',[0.5,1.5]);
    [probinfo,probcounter] = crashtestCheckOutput(vs,ls,exitflag,lams,lame,sprintf('\n### %d: dsTol: [0.5,1.5] ###\n',i),probinfo,probcounter);
    %% descale
    fprintf('\n### %d: scaling: staticdscale ###\n',i);
    [vs,ls,exitflag] = continuation(funJacoTest{i},v0,lams,lame,ds0,'dsMax',dsMax,'display','off','scaling','staticdscale');
    [probinfo,probcounter] = crashtestCheckOutput(vs,ls,exitflag,lams,lame,sprintf('\n### %d: scaling: staticdscale ###\n',i),probinfo,probcounter);
    fprintf('\n### %d: scaling: dynamicdscale ###\n',i);
    [vs,ls,exitflag] = continuation(funJacoTest{i},v0,lams,lame,ds0,'dsMax',dsMax,'display','off','scaling','dynamicdscale');
    [probinfo,probcounter] = crashtestCheckOutput(vs,ls,exitflag,lams,lame,sprintf('\n### %d: scaling: dynamicdscale ###\n',i),probinfo,probcounter);
    %% l0:
    fprintf('\n### %d: l0: (lEnd-lStart)/2+lStart ###\n',i);
    [vs,ls,exitflag] = continuation(funJacoTest{i},v0,lams,lame,ds0,'dsMax',dsMax,'display','off','l0',(lame-lams)/2+lams);
    [probinfo,probcounter] = crashtestCheckOutput(vs,ls,exitflag,lams,lame,sprintf('\n### %d: l0: (lEnd-lStart)/2+lStart ###\n',i),probinfo,probcounter);
    %% lTarget:
    fprintf('\n### %d: lTarget: (lEnd-lStart)/3+lStart ###\n',i);
    [vs,ls,exitflag] = continuation(funJacoTest{i},v0,lams,lame,ds0,'dsMax',dsMax,'dsMin',10^-6,'display','off','lTarget',(lame-lams)/3+lams);
    [probinfo,probcounter] = crashtestCheckOutput(vs,ls,exitflag,lams,lame,sprintf('\n### %d: lTarget: (lEnd-lStart)/3+lStart ###\n',i),probinfo,probcounter);
    %% lTarget:
    fprintf('\n### %d: alphaReverse: pi/4 ###\n',i);
    [vs,ls,exitflag] = continuation(funJacoTest{i},v0,lams,lame,ds0,'dsMax',dsMax,'display','off','alphaReverse',pi/4,'alphaReverseAutoMode','off');
    [probinfo,probcounter] = crashtestCheckOutput(vs,ls,exitflag,lams,lame,sprintf('\n### %d: alphaReverse: pi/4 ###\n',i),probinfo,probcounter);
    %% plot:
    fprintf('\n### %d: plot: on ###\n',i);
    [vs,ls,exitflag] = continuation(funJacoTest{i},v0,lams,lame,ds0,'dsMax',dsMax,'display','off','plot','on');
    [probinfo,probcounter] = crashtestCheckOutput(vs,ls,exitflag,lams,lame,sprintf('\n### %d: plot: on ###\n',i),probinfo,probcounter);
    fprintf('\n### %d: plot: on & bifurcation: mark ###\n',i);
    [vs,ls,exitflag] = continuation(funJacoTest{i},v0,lams,lame,ds0,'dsMax',dsMax,'display','off','plot','on','bifurcation','mark');
    [probinfo,probcounter] = crashtestCheckOutput(vs,ls,exitflag,lams,lame,sprintf('\n### %d: plot: on & bifurcation: mark ###\n',i),probinfo,probcounter);
    fprintf('\n### %d: plot: on & bifurcation: determine ###\n',i);
    [vs,ls,exitflag] = continuation(funJacoTest{i},v0,lams,lame,ds0,'dsMax',dsMax,'display','off','plot','on','bifurcation','determine');
    [probinfo,probcounter] = crashtestCheckOutput(vs,ls,exitflag,lams,lame,sprintf('\n### %d: plot: on & bifurcation: determine ###\n',i),probinfo,probcounter);
    fprintf('\n### %d: plot: on & bifurcation: trace ###\n',i);
    [vs,ls,exitflag] = continuation(funJacoTest{i},v0,lams,lame,ds0,'dsMax',dsMax,'display','off','plot','on','bifurcation','trace');
    [probinfo,probcounter] = crashtestCheckOutput(vs,ls,exitflag,lams,lame,sprintf('\n### %d: plot: on & bifurcation: trace ###\n',i),probinfo,probcounter);
    fprintf('\n### %d: plotVarsIndex: 1 ###\n',i);
    [vs,ls,exitflag] = continuation(funJacoTest{i},v0,lams,lame,ds0,'dsMax',dsMax,'display','off','plot','on','plotOptions',plot.PlotOptions('plotVarsIndex',1));
    [probinfo,probcounter] = crashtestCheckOutput(vs,ls,exitflag,lams,lame,sprintf('\n### %d: plotVarsIndex: 1 ###\n',i),probinfo,probcounter);
    %% includeReverse:
    fprintf('\n### %d: includeReverse: on ###\n',i);
    [vs,ls,exitflag] = continuation(funJacoTest{i},v0,lams,lame,ds0,'dsMax',dsMax,'display','off','includeReverse','on');
    [probinfo,probcounter] = crashtestCheckOutput(vs,ls,exitflag,lams,lame,sprintf('\n### %d: includeReverse: on ###\n',i),probinfo,probcounter);
    %% predictorTangential:
    fprintf('\n### %d: predictorTangential ###\n',i);
    [vs,ls,exitflag] = continuation(funJacoTest{i},v0,lams,lame,ds0,'dsMax',dsMax,'display','off','predictor','tangential');
    [probinfo,probcounter] = crashtestCheckOutput(vs,ls,exitflag,lams,lame,sprintf('\n### %d: predictorTangential ###\n',i),probinfo,probcounter);
    %% predictorPolynomialDegree:
    fprintf('\n### %d: predictorPolynomialDegree: 2 ###\n',i);
    [vs,ls,exitflag] = continuation(funJacoTest{i},v0,lams,lame,ds0,'dsMax',dsMax,'display','off','predictor','polynomial','predictorPolynomialDegree',2);
    [probinfo,probcounter] = crashtestCheckOutput(vs,ls,exitflag,lams,lame,sprintf('\n### %d: predictorPolynomialDegree: 2 ###\n',i),probinfo,probcounter);
    %% predictorPolynomialFit:
    fprintf('\n### %d: predictorPolynomialFit: 4 ###\n',i);
    [vs,ls,exitflag] = continuation(funJacoTest{i},v0,lams,lame,ds0,'dsMax',dsMax,'display','off','predictor','polynomial','predictorPolynomialFit',4);
    [probinfo,probcounter] = crashtestCheckOutput(vs,ls,exitflag,lams,lame,sprintf('\n### %d: predictorPolynomialFit: 4 ###\n',i),probinfo,probcounter);
    %% predictorPolynomialAdaptive:
    fprintf('\n### %d: predictorPolynomialAdaptive: nt = 4, nf = 5 ###\n',i);
    [vs,ls,exitflag] = continuation(funJacoTest{i},v0,lams,lame,ds0,'dsMax',dsMax,'display','off','predictor','polynomial','predictorPolynomialAdaptive','on','predictorPolynomialDegree',4,'predictorPolynomialFit',5);
    [probinfo,probcounter] = crashtestCheckOutput(vs,ls,exitflag,lams,lame,sprintf('\n### %d: predictorPolynomialAdaptive: nt = 4, nf = 5 ###\n',i),probinfo,probcounter);
    %% predictorPolynomialAdaptive:
    fprintf('\n### %d: predictorSolver: on ###\n',i);
    [vs,ls,exitflag] = continuation(funJacoTest{i},v0,lams,lame,ds0,'dsMax',dsMax,'display','off','predictorSolver','on');
    [probinfo,probcounter] = crashtestCheckOutput(vs,ls,exitflag,lams,lame,sprintf('\n### %d: predictorSolver: on ###\n',i),probinfo,probcounter);
    %% stepSizeControl:
    stepsizeMethods = ["angleChange", "angleCustom", "contraction" , "error", "errorAlt", "fayezioghani",...
        "iterationsExponential", "iterationsPolynomial", "multiplicative", "multiplicativeAlt", "pidCustom",...
        "pidValli", "szyszkowski", "yoon"];
    for stepsizeMethod = stepsizeMethods
        fprintf('\n### %d: stepSizeControl: %s ###\n',i, stepsizeMethod);
        [vs,ls,exitflag] = continuation(funJacoTest{i},v0,lams,lame,ds0,'dsMax',dsMax,'display','off','stepSizeControl',stepsizeMethod);
        [probinfo,probcounter] = crashtestCheckOutput(vs,ls,exitflag,lams,lame,sprintf('\n### %d: stepSizeControl: %s ###\n',i,stepsizeMethod),probinfo,probcounter);
    end
end
clc;
if isempty(probinfo)
    fprintf('###################################\n###################################\n### ### ends without errors ### ###\n###################################\n###################################\n');
else
    fprintf('###################################\n###################################\n### ### ends with problems  ### ###\n###################################\n###################################\n');
    fprintf('\n%s',probinfo,probcounter);
    fprintf('\n###################################\n###################################\n### ### ends with problems  ### ###\n###################################\n###################################\n');
    disp(['problems occured: ' num2str(probcounter)])
end
disp(['simulation time: ' num2str(toc)])