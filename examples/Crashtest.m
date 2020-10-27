clear all;
close all;
clc;
tic
addpath('..\src');
addpath('test_cases');
addpath('crashtest_aux');
probinfo = [];
probcounter = 0;
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
    probinfo = crashtest_check_output(vs,ls,exitflag,lams,lame,sprintf('### %d: standard-config ###\n',i),probinfo,probcounter); 
    %% display off:
    fprintf('\n### %d: display: off ###\n',i);
    [vs,ls,exitflag] = continuation(fun_jaco_test{i},v0,lams,lame,ds0,'display','off');
    [probinfo,probcounter] = crashtest_check_output(vs,ls,exitflag,lams,lame,sprintf('### %d: display: off ###\n',i),probinfo,probcounter); 
    %% homotopy:
    fprintf('\n### %d: homotopy: off ###\n',i);
    [vs,ls,exitflag] = continuation(fun_jaco_test{i},v0,lams,lame,ds0,'display','off','homotopy','off');
    [probinfo,probcounter] = crashtest_check_output(vs,ls,exitflag,lams,lame,sprintf('### %d: homotopy: off ###\n',i),probinfo,probcounter); 
    fprintf('\n### %d: homotopy: on ###\n',i);
    [vs,ls,exitflag] = continuation(fun_jaco_test{i},v0,lams,lame,ds0,'display','off','homotopy','on');
    [probinfo,probcounter] = crashtest_check_output(vs,ls,exitflag,lams,lame,sprintf('### %d: homotopy: on ###\n',i),probinfo,probcounter); 
    fprintf('\n### %d: homotopy: fix ###\n',i);
    [vs,ls,exitflag] = continuation(fun_jaco_test{i},v0,lams,lame,ds0,'display','off','homotopy','fix');
    [probinfo,probcounter] = crashtest_check_output(vs,ls,exitflag,lams,lame,sprintf('### %d: homotopy: fix ###\n',i),probinfo,probcounter); 
    fprintf('\n### %d: homotopy: newton ###\n',i);
    [vs,ls,exitflag] = continuation(fun_jaco_test{i},v0,lams,lame,ds0,'display','off','homotopy','newton');
    [probinfo,probcounter] = crashtest_check_output(vs,ls,exitflag,lams,lame,sprintf('### %d: homotopy: newton ###\n',i),probinfo,probcounter); 
    fprintf('\n### %d: homotopy: newton ###\n',i);
    [vs,ls,exitflag] = continuation(fun_jaco_test{i},v0,lams,lame,ds0,'display','off','homotopy','f2');
    [probinfo,probcounter] = crashtest_check_output(vs,ls,exitflag,lams,lame,sprintf('### %d: homotopy: newton ###\n',i),probinfo,probcounter); 
    %% solver:
    fprintf('\n### %d: solver: fsolve ###\n',i);
    [vs,ls,exitflag] = continuation(fun_jaco_test{i},v0,lams,lame,ds0,'display','off','solver','fsolve');
    [probinfo,probcounter] = crashtest_check_output(vs,ls,exitflag,lams,lame,sprintf('### %d: solver: fsolve ###\n',i),probinfo,probcounter); 
    fprintf('\n### %d: solver: fmincon ###\n',i);
    [vs,ls,exitflag] = continuation(fun_jaco_test{i},v0,lams,lame,ds0,'display','off','solver','fmincon');
    [probinfo,probcounter] = crashtest_check_output(vs,ls,exitflag,lams,lame,sprintf('### %d: solver: fmincon ###\n',i),probinfo,probcounter); 
    fprintf('\n### %d: solver: lsqnonlin ###\n',i);
    [vs,ls,exitflag] = continuation(fun_jaco_test{i},v0,lams,lame,ds0,'display','on','solver','lsqnonlin');
    [probinfo,probcounter] = crashtest_check_output(vs,ls,exitflag,lams,lame,sprintf('### %d: solver: lsqnonlin ###\n',i),probinfo,probcounter); 
    fprintf('\n### %d: solver: newton ###\n',i);
    [vs,ls,exitflag] = continuation(fun_jaco_test{i},v0,lams,lame,ds0,'display','off','solver','newton');
    [probinfo,probcounter] = crashtest_check_output(vs,ls,exitflag,lams,lame,sprintf('### %d: solver: newton ###\n',i),probinfo,probcounter); 
    %% arc-length:
    fprintf('\n### %d: arc-length: sphere ###\n',i);
    [vs,ls,exitflag] = continuation(fun_jaco_test{i},v0,lams,lame,ds0,'display','off','arclength','sphere');
    [probinfo,probcounter] = crashtest_check_output(vs,ls,exitflag,lams,lame,sprintf('### %d: arc-length: sphere ###\n',i),probinfo,probcounter); 
    fprintf('\n### %d: arc-length: linear ###\n',i);
    [vs,ls,exitflag] = continuation(fun_jaco_test{i},v0,lams,lame,ds0,'display','off','arclength','linear');
    [probinfo,probcounter] = crashtest_check_output(vs,ls,exitflag,lams,lame,sprintf('### %d: arc-length: linear ###\n',i),probinfo,probcounter); 
    fprintf('\n### %d: arc-length: ellipsoid ###\n',i);
    [vs,ls,exitflag] = continuation(fun_jaco_test{i},v0,lams,lame,ds0,'display','off','arclength','ellipsoid');
    [probinfo,probcounter] = crashtest_check_output(vs,ls,exitflag,lams,lame,sprintf('### %d: arc-length: ellipsoid ###\n',i),probinfo,probcounter); 
    fprintf('\n### %d: arc-length: ellipsoid2 ###\n',i);
    [vs,ls,exitflag] = continuation(fun_jaco_test{i},v0,lams,lame,ds0,'display','off','arclength','ellipsoid2');
    [probinfo,probcounter] = crashtest_check_output(vs,ls,exitflag,lams,lame,sprintf('### %d: arc-length: ellipsoid2 ###\n',i),probinfo,probcounter); 
    fprintf('\n### %d: arc-length: unique ###\n',i);
    [vs,ls,exitflag] = continuation(fun_jaco_test{i},v0,lams,lame,ds0,'display','off','arclength','unique');
    [probinfo,probcounter] = crashtest_check_output(vs,ls,exitflag,lams,lame,sprintf('### %d: arc-length: unique ###\n',i),probinfo,probcounter); 
    %% deflation:
    fprintf('\n### %d: deflation: on ###\n',i);
    [vs,ls,exitflag] = continuation(fun_jaco_test{i},v0,lams,lame,ds0,'display','off','deflation','on');
    [probinfo,probcounter] = crashtest_check_output(vs,ls,exitflag,lams,lame,sprintf('### %d: deflation: on ###\n',i),probinfo,probcounter); 
    fprintf('\n### %d: deflation: off ###\n',i);
    [vs,ls,exitflag] = continuation(fun_jaco_test{i},v0,lams,lame,ds0,'display','off','deflation','off');
    [probinfo,probcounter] = crashtest_check_output(vs,ls,exitflag,lams,lame,sprintf('### %d: deflation: off ###\n',i),probinfo,probcounter); 
    %% bifurcation:
    fprintf('\n### %d: bifurcation: mark ###\n',i);
    [vs,ls,exitflag] = continuation(fun_jaco_test{i},v0,lams,lame,ds0,'display','off','bifurcation','mark');
    [probinfo,probcounter] = crashtest_check_output(vs,ls,exitflag,lams,lame,sprintf('### %d: bifurcation: mark ###\n',i),probinfo,probcounter); 
    fprintf('\n### %d: bifurcation: determine ###\n',i);
    [vs,ls,exitflag] = continuation(fun_jaco_test{i},v0,lams,lame,ds0,'display','off','bifurcation','determine');
    [probinfo,probcounter] = crashtest_check_output(vs,ls,exitflag,lams,lame,sprintf('### %d: bifurcation: determine ###\n',i),probinfo,probcounter); 
    fprintf('\n### %d: bifurcation: trace ###\n',i);
    [vs,ls,exitflag] = continuation(fun_jaco_test{i},v0,lams,lame,ds0,'display','off','bifurcation','trace');
    [probinfo,probcounter] = crashtest_check_output(vs,ls,exitflag,lams,lame,sprintf('### %d: bifurcation: trace ###\n',i),probinfo,probcounter); 
    %% n_iter_opt:
    fprintf('\n### %d: n_iter_opt: 5 ###\n',i);
    [vs,ls,exitflag] = continuation(fun_jaco_test{i},v0,lams,lame,ds0,'display','off','n_iter_opt',5);
    [probinfo,probcounter] = crashtest_check_output(vs,ls,exitflag,lams,lame,sprintf('### %d: n_iter_opt: 5 ###\n',i),probinfo,probcounter); 
    %% ds_max:
    fprintf('\n### %d: ds_max: 2*ds0 ###\n',i);
    [vs,ls,exitflag] = continuation(fun_jaco_test{i},v0,lams,lame,ds0,'display','off','ds_max',2*ds0);
    [probinfo,probcounter] = crashtest_check_output(vs,ls,exitflag,lams,lame,sprintf('### %d: ds_max: 2*ds0 ###\n',i),probinfo,probcounter); 
    %% l_0:
    fprintf('\n### %d: l_0: (l_end-l_start)/2+l_start ###\n',i);
    [vs,ls,exitflag] = continuation(fun_jaco_test{i},v0,lams,lame,ds0,'display','off','l_0',(lame-lams)/2+lams);
    [probinfo,probcounter] = crashtest_check_output(vs,ls,exitflag,lams,lame,sprintf('### %d: l_0: (l_end-l_start)/2+l_start ###\n',i),probinfo,probcounter); 
    %% l_target:
    fprintf('\n### %d: l_target: (l_end-l_start)/2+l_start ###\n',i);
    [vs,ls,exitflag] = continuation(fun_jaco_test{i},v0,lams,lame,ds0,'display','off','l_target',(lame-lams)/2+lams);
    [probinfo,probcounter] = crashtest_check_output(vs,ls,exitflag,lams,lame,sprintf('### %d: l_target: (l_end-l_start)/2+l_start ###\n',i),probinfo,probcounter); 
    %% l_target:
    fprintf('\n### %d: alpha_reverse: pi/4 ###\n',i);
    [vs,ls,exitflag] = continuation(fun_jaco_test{i},v0,lams,lame,ds0,'display','off','alpha_reverse',pi/4);
    [probinfo,probcounter] = crashtest_check_output(vs,ls,exitflag,lams,lame,sprintf('### %d: alpha_reverse: pi/4 ###\n',i),probinfo,probcounter); 
    %% unique:
    fprintf('\n### %d: unique: on ###\n',i);
    [vs,ls,exitflag] = continuation(fun_jaco_test{i},v0,lams,lame,ds0,'display','off','unique','on');
    [probinfo,probcounter] = crashtest_check_output(vs,ls,exitflag,lams,lame,sprintf('### %d: unique: on ###\n',i),probinfo,probcounter); 
    %% plot:
    fprintf('\n### %d: plot: on ###\n',i);
    [vs,ls,exitflag] = continuation(fun_jaco_test{i},v0,lams,lame,ds0,'display','off','plot','on');
    [probinfo,probcounter] = crashtest_check_output(vs,ls,exitflag,lams,lame,sprintf('### %d: plot: on ###\n',i),probinfo,probcounter); 
    fprintf('\n### %d: plot_vars_index: 1 ###\n',i);
    [vs,ls,exitflag] = continuation(fun_jaco_test{i},v0,lams,lame,ds0,'display','off','plot','on','plot_vars_index',1);
    [probinfo,probcounter] = crashtest_check_output(vs,ls,exitflag,lams,lame,sprintf('### %d: plot_vars_index: 1 ###\n',i),probinfo,probcounter); 
    %% include_reverse:
    fprintf('\n### %d: include_reverse: on ###\n',i);
    [vs,ls,exitflag] = continuation(fun_jaco_test{i},v0,lams,lame,ds0,'display','off','include_reverse','on');
    [probinfo,probcounter] = crashtest_check_output(vs,ls,exitflag,lams,lame,sprintf('### %d: include_reverse: on ###\n',i),probinfo,probcounter); 
    %% predictor_taylor:
    fprintf('\n### %d: predictor_taylor: 2 ###\n',i);
    [vs,ls,exitflag] = continuation(fun_jaco_test{i},v0,lams,lame,ds0,'display','off','predictor_taylor',2);
    [probinfo,probcounter] = crashtest_check_output(vs,ls,exitflag,lams,lame,sprintf('### %d: predictor_taylor: 2 ###\n',i),probinfo,probcounter); 
    %% predictor_fit:
    fprintf('\n### %d: predictor_fit: 4 ###\n',i);
    [vs,ls,exitflag] = continuation(fun_jaco_test{i},v0,lams,lame,ds0,'display','off','predictor_fit',4);
    [probinfo,probcounter] = crashtest_check_output(vs,ls,exitflag,lams,lame,sprintf('### %d: predictor_fit: 4 ###\n',i),probinfo,probcounter); 
    %% predictor_adaptive:
    fprintf('\n### %d: predictor_adaptive: nt = 4, nf = 5 ###\n',i);
    [vs,ls,exitflag] = continuation(fun_jaco_test{i},v0,lams,lame,ds0,'display','off','predictor_adaptive','on','predictor_taylor',4,'predictor_fit',5);
    [probinfo,probcounter] = crashtest_check_output(vs,ls,exitflag,lams,lame,sprintf('### %d: predictor_adaptive: nt = 4, nf = 5 ###\n',i),probinfo,probcounter); 
    %% step_size_control:
    fprintf('\n### %d: step_size_control: angle ###\n',i);
    [vs,ls,exitflag] = continuation(fun_jaco_test{i},v0,lams,lame,ds0,'display','off','step_size_control','angle');
    [probinfo,probcounter] = crashtest_check_output(vs,ls,exitflag,lams,lame,sprintf('### %d: step_size_control: angle ###\n',i),probinfo,probcounter); 
    fprintf('\n### %d: step_size_control: pid ###\n',i);
    [vs,ls,exitflag] = continuation(fun_jaco_test{i},v0,lams,lame,ds0,'display','off','step_size_control','pid');
    [probinfo,probcounter] = crashtest_check_output(vs,ls,exitflag,lams,lame,sprintf('### %d: step_size_control: pid ###\n',i),probinfo,probcounter); 
    fprintf('\n### %d: step_size_control: curvature ###\n',i);
    [vs,ls,exitflag] = continuation(fun_jaco_test{i},v0,lams,lame,ds0,'display','off','step_size_control','curvature');
    [probinfo,probcounter] = crashtest_check_output(vs,ls,exitflag,lams,lame,sprintf('### %d: step_size_control: curvature ###\n',i),probinfo,probcounter); 
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