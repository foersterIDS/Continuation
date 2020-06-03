clear all;
clc;
tic
addpath('..\src');
addpath('test_cases');
%% Testfunktion #04
% Point of intersection of circle and sin(radius)-scaled exponential function with radius as parameter
% testfun04;
%% Testfunktion #05
% function with bifurcation
testfun05;
%% Crashtest:
for i=1:2
    % i=1: with jacobian, i=2: without jacobian
    %% standard-config.:
    fprintf('\n### %d: standard-config ###\n',i);
    [vs,ls,exitflag] = continuation(fun_jaco_test{i},v0,lams,lame,ds0);
    %% display off:
    fprintf('\n### %d: display: off ###\n',i);
    [vs,ls,exitflag] = continuation(fun_jaco_test{i},v0,lams,lame,ds0,'display','off');
    %% homotopy:
    fprintf('\n### %d: homotopy: off ###\n',i);
    [vs,ls,exitflag] = continuation(fun_jaco_test{i},v0,lams,lame,ds0,'display','off','homotopy','off');
    fprintf('\n### %d: homotopy: on ###\n',i);
    [vs,ls,exitflag] = continuation(fun_jaco_test{i},v0,lams,lame,ds0,'display','off','homotopy','on');
    fprintf('\n### %d: homotopy: fix ###\n',i);
    [vs,ls,exitflag] = continuation(fun_jaco_test{i},v0,lams,lame,ds0,'display','off','homotopy','fix');
    fprintf('\n### %d: homotopy: newton ###\n',i);
    [vs,ls,exitflag] = continuation(fun_jaco_test{i},v0,lams,lame,ds0,'display','off','homotopy','newton');
    %% solver:
    fprintf('\n### %d: solver: fsolve ###\n',i);
    [vs,ls,exitflag] = continuation(fun_jaco_test{i},v0,lams,lame,ds0,'display','off','solver','fsolve');
    fprintf('\n### %d: solver: fmincon ###\n',i);
    [vs,ls,exitflag] = continuation(fun_jaco_test{i},v0,lams,lame,ds0,'display','off','solver','fmincon');
    fprintf('\n### %d: solver: lsqnonlin ###\n',i);
    [vs,ls,exitflag] = continuation(fun_jaco_test{i},v0,lams,lame,ds0,'display','off','solver','lsqnonlin');
    fprintf('\n### %d: solver: newton ###\n',i);
    [vs,ls,exitflag] = continuation(fun_jaco_test{i},v0,lams,lame,ds0,'display','off','solver','newton');
    %% arc-length:
    fprintf('\n### %d: arc-length: sphere ###\n',i);
    [vs,ls,exitflag] = continuation(fun_jaco_test{i},v0,lams,lame,ds0,'display','off','arclength','sphere');
    fprintf('\n### %d: arc-length: linear ###\n',i);
    [vs,ls,exitflag] = continuation(fun_jaco_test{i},v0,lams,lame,ds0,'display','off','arclength','linear');
    fprintf('\n### %d: arc-length: ellipsoid ###\n',i);
    [vs,ls,exitflag] = continuation(fun_jaco_test{i},v0,lams,lame,ds0,'display','off','arclength','ellipsoid');
    %% deflation:
    fprintf('\n### %d: deflation: on ###\n',i);
    [vs,ls,exitflag] = continuation(fun_jaco_test{i},v0,lams,lame,ds0,'display','off','deflation','on');
    fprintf('\n### %d: deflation: off ###\n',i);
    [vs,ls,exitflag] = continuation(fun_jaco_test{i},v0,lams,lame,ds0,'display','off','deflation','off');
    %% bifurcation:
    fprintf('\n### %d: bifurcation: mark ###\n',i);
    [vs,ls,exitflag] = continuation(fun_jaco_test{i},v0,lams,lame,ds0,'display','off','bifurcation','mark');
    fprintf('\n### %d: bifurcation: determine ###\n',i);
    [vs,ls,exitflag] = continuation(fun_jaco_test{i},v0,lams,lame,ds0,'display','off','bifurcation','determine');
    fprintf('\n### %d: bifurcation: trace ###\n',i);
    [vs,ls,exitflag] = continuation(fun_jaco_test{i},v0,lams,lame,ds0,'display','off','bifurcation','trace');
end
fprintf('\n########################\n########################\n### ### success! ### ###\n########################\n########################\n');
disp(['simulation time: ' num2str(toc)])