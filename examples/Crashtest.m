clear all;
clc;
tic
addpath('..\src');
addpath('test_cases');
addpath('crashtest_aux');
probinfo = [];
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
    probinfo = crashtest_check_output(vs,ls,exitflag,lams,lame,sprintf('### %d: standard-config ###\n',i),probinfo);
    %% display off:
    fprintf('\n### %d: display: off ###\n',i);
    [vs,ls,exitflag] = continuation(fun_jaco_test{i},v0,lams,lame,ds0,'display','off');
    probinfo = crashtest_check_output(vs,ls,exitflag,lams,lame,sprintf('### %d: display: off ###\n',i),probinfo);
    %% homotopy:
    fprintf('\n### %d: homotopy: off ###\n',i);
    [vs,ls,exitflag] = continuation(fun_jaco_test{i},v0,lams,lame,ds0,'display','off','homotopy','off');
    probinfo = crashtest_check_output(vs,ls,exitflag,lams,lame,sprintf('### %d: homotopy: off ###\n',i),probinfo);
    fprintf('\n### %d: homotopy: on ###\n',i);
    [vs,ls,exitflag] = continuation(fun_jaco_test{i},v0,lams,lame,ds0,'display','off','homotopy','on');
    probinfo = crashtest_check_output(vs,ls,exitflag,lams,lame,sprintf('### %d: homotopy: on ###\n',i),probinfo);
    fprintf('\n### %d: homotopy: fix ###\n',i);
    [vs,ls,exitflag] = continuation(fun_jaco_test{i},v0,lams,lame,ds0,'display','off','homotopy','fix');
    probinfo = crashtest_check_output(vs,ls,exitflag,lams,lame,sprintf('### %d: homotopy: fix ###\n',i),probinfo);
    fprintf('\n### %d: homotopy: newton ###\n',i);
    [vs,ls,exitflag] = continuation(fun_jaco_test{i},v0,lams,lame,ds0,'display','off','homotopy','newton');
    probinfo = crashtest_check_output(vs,ls,exitflag,lams,lame,sprintf('### %d: homotopy: newton ###\n',i),probinfo);
    fprintf('\n### %d: homotopy: newton ###\n',i);
    [vs,ls,exitflag] = continuation(fun_jaco_test{i},v0,lams,lame,ds0,'display','off','homotopy','f2');
    probinfo = crashtest_check_output(vs,ls,exitflag,lams,lame,sprintf('### %d: homotopy: newton ###\n',i),probinfo);
    %% solver:
    fprintf('\n### %d: solver: fsolve ###\n',i);
    [vs,ls,exitflag] = continuation(fun_jaco_test{i},v0,lams,lame,ds0,'display','off','solver','fsolve');
    probinfo = crashtest_check_output(vs,ls,exitflag,lams,lame,sprintf('### %d: solver: fsolve ###\n',i),probinfo);
    fprintf('\n### %d: solver: fmincon ###\n',i);
    [vs,ls,exitflag] = continuation(fun_jaco_test{i},v0,lams,lame,ds0,'display','off','solver','fmincon');
    probinfo = crashtest_check_output(vs,ls,exitflag,lams,lame,sprintf('### %d: solver: fmincon ###\n',i),probinfo);
    fprintf('\n### %d: solver: lsqnonlin ###\n',i);
    [vs,ls,exitflag] = continuation(fun_jaco_test{i},v0,lams,lame,ds0,'display','on','solver','lsqnonlin');
    probinfo = crashtest_check_output(vs,ls,exitflag,lams,lame,sprintf('### %d: solver: lsqnonlin ###\n',i),probinfo);
    fprintf('\n### %d: solver: newton ###\n',i);
    [vs,ls,exitflag] = continuation(fun_jaco_test{i},v0,lams,lame,ds0,'display','off','solver','newton');
    probinfo = crashtest_check_output(vs,ls,exitflag,lams,lame,sprintf('### %d: solver: newton ###\n',i),probinfo);
    %% arc-length:
    fprintf('\n### %d: arc-length: sphere ###\n',i);
    [vs,ls,exitflag] = continuation(fun_jaco_test{i},v0,lams,lame,ds0,'display','off','arclength','sphere');
    probinfo = crashtest_check_output(vs,ls,exitflag,lams,lame,sprintf('### %d: arc-length: sphere ###\n',i),probinfo);
    fprintf('\n### %d: arc-length: linear ###\n',i);
    [vs,ls,exitflag] = continuation(fun_jaco_test{i},v0,lams,lame,ds0,'display','off','arclength','linear');
    probinfo = crashtest_check_output(vs,ls,exitflag,lams,lame,sprintf('### %d: arc-length: linear ###\n',i),probinfo);
    fprintf('\n### %d: arc-length: ellipsoid ###\n',i);
    [vs,ls,exitflag] = continuation(fun_jaco_test{i},v0,lams,lame,ds0,'display','off','arclength','ellipsoid');
    probinfo = crashtest_check_output(vs,ls,exitflag,lams,lame,sprintf('### %d: arc-length: ellipsoid ###\n',i),probinfo);
    fprintf('\n### %d: arc-length: ellipsoid2 ###\n',i);
    [vs,ls,exitflag] = continuation(fun_jaco_test{i},v0,lams,lame,ds0,'display','off','arclength','ellipsoid2');
    probinfo = crashtest_check_output(vs,ls,exitflag,lams,lame,sprintf('### %d: arc-length: ellipsoid2 ###\n',i),probinfo);
    fprintf('\n### %d: arc-length: unique ###\n',i);
    [vs,ls,exitflag] = continuation(fun_jaco_test{i},v0,lams,lame,ds0,'display','off','arclength','unique');
    probinfo = crashtest_check_output(vs,ls,exitflag,lams,lame,sprintf('### %d: arc-length: unique ###\n',i),probinfo);
    %% deflation:
    fprintf('\n### %d: deflation: on ###\n',i);
    [vs,ls,exitflag] = continuation(fun_jaco_test{i},v0,lams,lame,ds0,'display','off','deflation','on');
    probinfo = crashtest_check_output(vs,ls,exitflag,lams,lame,sprintf('### %d: deflation: on ###\n',i),probinfo);
    fprintf('\n### %d: deflation: off ###\n',i);
    [vs,ls,exitflag] = continuation(fun_jaco_test{i},v0,lams,lame,ds0,'display','off','deflation','off');
    probinfo = crashtest_check_output(vs,ls,exitflag,lams,lame,sprintf('### %d: deflation: off ###\n',i),probinfo);
    %% bifurcation:
    fprintf('\n### %d: bifurcation: mark ###\n',i);
    [vs,ls,exitflag] = continuation(fun_jaco_test{i},v0,lams,lame,ds0,'display','off','bifurcation','mark');
    probinfo = crashtest_check_output(vs,ls,exitflag,lams,lame,sprintf('### %d: bifurcation: mark ###\n',i),probinfo);
    fprintf('\n### %d: bifurcation: determine ###\n',i);
    [vs,ls,exitflag] = continuation(fun_jaco_test{i},v0,lams,lame,ds0,'display','off','bifurcation','determine');
    probinfo = crashtest_check_output(vs,ls,exitflag,lams,lame,sprintf('### %d: bifurcation: determine ###\n',i),probinfo);
    fprintf('\n### %d: bifurcation: trace ###\n',i);
    [vs,ls,exitflag] = continuation(fun_jaco_test{i},v0,lams,lame,ds0,'display','off','bifurcation','trace');
    probinfo = crashtest_check_output(vs,ls,exitflag,lams,lame,sprintf('### %d: bifurcation: trace ###\n',i),probinfo);
    %% n_iter_opt:
    fprintf('\n### %d: n_iter_opt: 5 ###\n',i);
    [vs,ls,exitflag] = continuation(fun_jaco_test{i},v0,lams,lame,ds0,'display','off','n_iter_opt',5);
    probinfo = crashtest_check_output(vs,ls,exitflag,lams,lame,sprintf('### %d: n_iter_opt: 5 ###\n',i),probinfo);
    %% ds_max:
    fprintf('\n### %d: ds_max: 2*ds0 ###\n',i);
    [vs,ls,exitflag] = continuation(fun_jaco_test{i},v0,lams,lame,ds0,'display','off','ds_max',2*ds0);
    probinfo = crashtest_check_output(vs,ls,exitflag,lams,lame,sprintf('### %d: ds_max: 2*ds0 ###\n',i),probinfo);
    %% l_0:
    fprintf('\n### %d: l_0: (l_end-l_start)/2+l_start ###\n',i);
    [vs,ls,exitflag] = continuation(fun_jaco_test{i},v0,lams,lame,ds0,'display','off','l_0',(lame-lams)/2+lams);
    probinfo = crashtest_check_output(vs,ls,exitflag,lams,lame,sprintf('### %d: l_0: (l_end-l_start)/2+l_start ###\n',i),probinfo);
    %% l_target:
    fprintf('\n### %d: l_target: (l_end-l_start)/2+l_start ###\n',i);
    [vs,ls,exitflag] = continuation(fun_jaco_test{i},v0,lams,lame,ds0,'display','off','l_target',(lame-lams)/2+lams);
    probinfo = crashtest_check_output(vs,ls,exitflag,lams,lame,sprintf('### %d: l_target: (l_end-l_start)/2+l_start ###\n',i),probinfo);
    %% l_target:
    fprintf('\n### %d: alpha_reverse: pi/4 ###\n',i);
    [vs,ls,exitflag] = continuation(fun_jaco_test{i},v0,lams,lame,ds0,'display','off','alpha_reverse',pi/4);
    probinfo = crashtest_check_output(vs,ls,exitflag,lams,lame,sprintf('### %d: alpha_reverse: pi/4 ###\n',i),probinfo);
    %% unique:
    fprintf('\n### %d: unique: on ###\n',i);
    [vs,ls,exitflag] = continuation(fun_jaco_test{i},v0,lams,lame,ds0,'display','off','unique','on');
    probinfo = crashtest_check_output(vs,ls,exitflag,lams,lame,sprintf('### %d: unique: on ###\n',i),probinfo);
    %% plot:
    fprintf('\n### %d: plot: on ###\n',i);
    [vs,ls,exitflag] = continuation(fun_jaco_test{i},v0,lams,lame,ds0,'display','off','plot','on');
    probinfo = crashtest_check_output(vs,ls,exitflag,lams,lame,sprintf('### %d: plot: on ###\n',i),probinfo);
    %% include_reverse:
    fprintf('\n### %d: include_reverse: on ###\n',i);
    [vs,ls,exitflag] = continuation(fun_jaco_test{i},v0,lams,lame,ds0,'display','off','include_reverse','on');
    probinfo = crashtest_check_output(vs,ls,exitflag,lams,lame,sprintf('### %d: include_reverse: on ###\n',i),probinfo);
    %% predictor_taylor:
    fprintf('\n### %d: predictor_taylor: 2 ###\n',i);
    [vs,ls,exitflag] = continuation(fun_jaco_test{i},v0,lams,lame,ds0,'display','off','predictor_taylor',2);
    probinfo = crashtest_check_output(vs,ls,exitflag,lams,lame,sprintf('### %d: predictor_taylor: 2 ###\n',i),probinfo);
    %% predictor_fit:
    fprintf('\n### %d: predictor_fit: 4 ###\n',i);
    [vs,ls,exitflag] = continuation(fun_jaco_test{i},v0,lams,lame,ds0,'display','off','predictor_fit',4);
    probinfo = crashtest_check_output(vs,ls,exitflag,lams,lame,sprintf('### %d: predictor_fit: 4 ###\n',i),probinfo);
    %% predictor_adaptive:
    fprintf('\n### %d: predictor_adaptive: nt = 4, nf = 5 ###\n',i);
    [vs,ls,exitflag] = continuation(fun_jaco_test{i},v0,lams,lame,ds0,'display','off','predictor_adaptive','on','predictor_taylor',4,'predictor_fit',5);
    probinfo = crashtest_check_output(vs,ls,exitflag,lams,lame,sprintf('### %d: predictor_adaptive: nt = 4, nf = 5 ###\n',i),probinfo);
    %% step_size_control:
    fprintf('\n### %d: step_size_control: angle ###\n',i);
    [vs,ls,exitflag] = continuation(fun_jaco_test{i},v0,lams,lame,ds0,'display','off','step_size_control','angle');
    probinfo = crashtest_check_output(vs,ls,exitflag,lams,lame,sprintf('### %d: step_size_control: angle ###\n',i),probinfo);
    fprintf('\n### %d: step_size_control: pid ###\n',i);
    [vs,ls,exitflag] = continuation(fun_jaco_test{i},v0,lams,lame,ds0,'display','off','step_size_control','pid');
    probinfo = crashtest_check_output(vs,ls,exitflag,lams,lame,sprintf('### %d: step_size_control: pid ###\n',i),probinfo);
    fprintf('\n### %d: step_size_control: curvature ###\n',i);
    [vs,ls,exitflag] = continuation(fun_jaco_test{i},v0,lams,lame,ds0,'display','off','step_size_control','curvature');
    probinfo = crashtest_check_output(vs,ls,exitflag,lams,lame,sprintf('### %d: step_size_control: curvature ###\n',i),probinfo);
end
clc;
if isempty(probinfo)
    fprintf('###################################\n###################################\n### ### ends without errors ### ###\n###################################\n###################################\n');
else
    fprintf('###################################\n###################################\n### ### ends with problems  ### ###\n###################################\n###################################\n');
    fprintf('\n%s',probinfo);
    fprintf('\n###################################\n###################################\n### ### ends with problems  ### ###\n###################################\n###################################\n');
end
disp(['simulation time: ' num2str(toc)])