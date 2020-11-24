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
testfun08;
%% Crashtest:
for i=1:2
    % i=1: with jacobian, i=2: without jacobian
    %% standard-config.:
    fprintf('\n### %d: standard-config ###\n',i);
    [vs,ls,exitflag] = continuation(fun_jaco_test{i},v0,lams,lame,ds0,'ds_max',ds_max);
    probinfo = crashtest_check_output(vs,ls,exitflag,lams,lame,sprintf('### %d: standard-config ###\n',i),probinfo,probcounter);
    %% display off:
    fprintf('\n### %d: display: off ###\n',i);
    [vs,ls,exitflag] = continuation(fun_jaco_test{i},v0,lams,lame,ds0,'ds_max',ds_max,'display','off');
    [probinfo,probcounter] = crashtest_check_output(vs,ls,exitflag,lams,lame,sprintf('### %d: display: off ###\n',i),probinfo,probcounter);
    %% homotopy:
    fprintf('\n### %d: homotopy: off ###\n',i);
    [vs,ls,exitflag] = continuation(fun_jaco_test{i},v0,lams,lame,ds0,'ds_max',ds_max,'display','off','homotopy','off');
    [probinfo,probcounter] = crashtest_check_output(vs,ls,exitflag,lams,lame,sprintf('### %d: homotopy: off ###\n',i),probinfo,probcounter);
    fprintf('\n### %d: homotopy: on ###\n',i);
    [vs,ls,exitflag] = continuation(fun_jaco_test{i},v0,lams,lame,ds0,'ds_max',ds_max,'display','off','homotopy','on');
    [probinfo,probcounter] = crashtest_check_output(vs,ls,exitflag,lams,lame,sprintf('### %d: homotopy: on ###\n',i),probinfo,probcounter);
    fprintf('\n### %d: homotopy: fix ###\n',i);
    [vs,ls,exitflag] = continuation(fun_jaco_test{i},v0,lams,lame,ds0,'ds_max',ds_max,'display','off','homotopy','fix');
    [probinfo,probcounter] = crashtest_check_output(vs,ls,exitflag,lams,lame,sprintf('### %d: homotopy: fix ###\n',i),probinfo,probcounter);
    fprintf('\n### %d: homotopy: newton ###\n',i);
    [vs,ls,exitflag] = continuation(fun_jaco_test{i},v0,lams,lame,ds0,'ds_max',ds_max,'display','off','homotopy','newton');
    [probinfo,probcounter] = crashtest_check_output(vs,ls,exitflag,lams,lame,sprintf('### %d: homotopy: newton ###\n',i),probinfo,probcounter);
    fprintf('\n### %d: homotopy: newton ###\n',i);
    [vs,ls,exitflag] = continuation(fun_jaco_test{i},v0,lams,lame,ds0,'ds_max',ds_max,'display','off','homotopy','f2');
    [probinfo,probcounter] = crashtest_check_output(vs,ls,exitflag,lams,lame,sprintf('### %d: homotopy: newton ###\n',i),probinfo,probcounter);
    %% solver:
    fprintf('\n### %d: solver: fsolve ###\n',i);
    [vs,ls,exitflag] = continuation(fun_jaco_test{i},v0,lams,lame,ds0,'ds_max',ds_max,'display','off','solver','fsolve');
    [probinfo,probcounter] = crashtest_check_output(vs,ls,exitflag,lams,lame,sprintf('### %d: solver: fsolve ###\n',i),probinfo,probcounter);
    fprintf('\n### %d: solver: lsqnonlin ###\n',i);
    [vs,ls,exitflag] = continuation(fun_jaco_test{i},v0,lams,lame,ds0,'ds_max',ds_max,'display','on','solver','lsqnonlin');
    [probinfo,probcounter] = crashtest_check_output(vs,ls,exitflag,lams,lame,sprintf('### %d: solver: lsqnonlin ###\n',i),probinfo,probcounter);
    fprintf('\n### %d: solver: newton ###\n',i);
    [vs,ls,exitflag] = continuation(fun_jaco_test{i},v0,lams,lame,ds0,'ds_max',ds_max,'display','off','solver','newton');
    [probinfo,probcounter] = crashtest_check_output(vs,ls,exitflag,lams,lame,sprintf('### %d: solver: newton ###\n',i),probinfo,probcounter);
    %% arc-length:
    fprintf('\n### %d: arc-length: sphere ###\n',i);
    [vs,ls,exitflag] = continuation(fun_jaco_test{i},v0,lams,lame,ds0,'ds_max',ds_max,'display','off','arclength','sphere');
    [probinfo,probcounter] = crashtest_check_output(vs,ls,exitflag,lams,lame,sprintf('### %d: arc-length: sphere ###\n',i),probinfo,probcounter);
    fprintf('\n### %d: arc-length: linear ###\n',i);
    [vs,ls,exitflag] = continuation(fun_jaco_test{i},v0,lams,lame,ds0,'ds_max',ds_max,'display','off','arclength','linear');
    [probinfo,probcounter] = crashtest_check_output(vs,ls,exitflag,lams,lame,sprintf('### %d: arc-length: linear ###\n',i),probinfo,probcounter);
    fprintf('\n### %d: arc-length: ellipsoid ###\n',i);
    [vs,ls,exitflag] = continuation(fun_jaco_test{i},v0,lams,lame,ds0,'ds_max',ds_max,'display','on','arclength','ellipsoid','predictor','tangential');
    [probinfo,probcounter] = crashtest_check_output(vs,ls,exitflag,lams,lame,sprintf('### %d: arc-length: ellipsoid ###\n',i),probinfo,probcounter);
    fprintf('\n### %d: arc-length: ellipsoid2 ###\n',i);
    [vs,ls,exitflag] = continuation(fun_jaco_test{i},v0,lams,lame,ds0,'ds_max',ds_max,'display','off','arclength','ellipsoid2');
    [probinfo,probcounter] = crashtest_check_output(vs,ls,exitflag,lams,lame,sprintf('### %d: arc-length: ellipsoid2 ###\n',i),probinfo,probcounter);
    fprintf('\n### %d: arc-length: unique ###\n',i);
    [vs,ls,exitflag] = continuation(fun_jaco_test{i},v0,lams,lame,ds0,'ds_max',ds_max,'display','off','arclength','unique');
    [probinfo,probcounter] = crashtest_check_output(vs,ls,exitflag,lams,lame,sprintf('### %d: arc-length: unique ###\n',i),probinfo,probcounter);
    %% deflation:
    fprintf('\n### %d: deflation: on ###\n',i);
    [vs,ls,exitflag] = continuation(fun_jaco_test{i},v0,lams,lame,ds0,'ds_max',ds_max,'display','off','deflation','on');
    [probinfo,probcounter] = crashtest_check_output(vs,ls,exitflag,lams,lame,sprintf('### %d: deflation: on ###\n',i),probinfo,probcounter);
    fprintf('\n### %d: deflation: off ###\n',i);
    [vs,ls,exitflag] = continuation(fun_jaco_test{i},v0,lams,lame,ds0,'ds_max',ds_max,'display','off','deflation','off');
    [probinfo,probcounter] = crashtest_check_output(vs,ls,exitflag,lams,lame,sprintf('### %d: deflation: off ###\n',i),probinfo,probcounter);
    %% bifurcation:
    fprintf('\n### %d: bifurcation: mark ###\n',i);
    [vs,ls,exitflag] = continuation(fun_jaco_test{i},v0,lams,lame,ds0,'ds_max',ds_max,'display','off','bifurcation','mark');
    [probinfo,probcounter] = crashtest_check_output(vs,ls,exitflag,lams,lame,sprintf('### %d: bifurcation: mark ###\n',i),probinfo,probcounter);
    fprintf('\n### %d: bifurcation: determine ###\n',i);
    [vs,ls,exitflag] = continuation(fun_jaco_test{i},v0,lams,lame,ds0,'ds_max',ds_max,'display','off','bifurcation','determine');
    [probinfo,probcounter] = crashtest_check_output(vs,ls,exitflag,lams,lame,sprintf('### %d: bifurcation: determine ###\n',i),probinfo,probcounter);
    fprintf('\n### %d: bifurcation: trace ###\n',i);
    [vs,ls,exitflag] = continuation(fun_jaco_test{i},v0,lams,lame,ds0,'ds_max',ds_max,'display','off','bifurcation','trace');
    [probinfo,probcounter] = crashtest_check_output(vs,ls,exitflag,lams,lame,sprintf('### %d: bifurcation: trace ###\n',i),probinfo,probcounter);
    fprintf('\n### %d: stop_on_bifurcation: on ###\n',i);
    [vs,ls,exitflag] = continuation(fun_jaco_test{i},v0,lams,lame,ds0,'ds_max',ds_max,'display','off','stop_on_bifurcation','on');
    [probinfo,probcounter] = crashtest_check_output(vs,ls,exitflag,lams,lame,sprintf('### %d: stop_on_bifurcation: on ###\n',i),probinfo,probcounter);
    %% n_iter_opt:
    fprintf('\n### %d: n_iter_opt: 5 ###\n',i);
    [vs,ls,exitflag] = continuation(fun_jaco_test{i},v0,lams,lame,ds0,'ds_max',ds_max,'display','off','n_iter_opt',5);
    [probinfo,probcounter] = crashtest_check_output(vs,ls,exitflag,lams,lame,sprintf('### %d: n_iter_opt: 5 ###\n',i),probinfo,probcounter);
    %% ds_min:
    fprintf('\n### %d: ds_min: ds0/2 ###\n',i);
    [vs,ls,exitflag] = continuation(fun_jaco_test{i},v0,lams,lame,ds0,'ds_max',ds_max,'display','off','ds_min',ds0/2);
    [probinfo,probcounter] = crashtest_check_output(vs,ls,exitflag,lams,lame,sprintf('### %d: ds_min: ds0/2 ###\n',i),probinfo,probcounter);
    %% descale
    fprintf('\n### %d: scaling: staticdscale ###\n',i);
    [vs,ls,exitflag] = continuation(fun_jaco_test{i},v0,lams,lame,ds0,'ds_max',ds_max,'display','off','scaling','staticdscale');
    [probinfo,probcounter] = crashtest_check_output(vs,ls,exitflag,lams,lame,sprintf('### %d: scaling: staticdscale ###\n',i),probinfo,probcounter);
    fprintf('\n### %d: scaling: dynamicdscale ###\n',i);
    [vs,ls,exitflag] = continuation(fun_jaco_test{i},v0,lams,lame,ds0,'ds_max',ds_max,'display','off','scaling','dynamicdscale');
    [probinfo,probcounter] = crashtest_check_output(vs,ls,exitflag,lams,lame,sprintf('### %d: scaling: dynamicdscale ###\n',i),probinfo,probcounter);
    %% l_0:
    fprintf('\n### %d: l_0: (l_end-l_start)/2+l_start ###\n',i);
    [vs,ls,exitflag] = continuation(fun_jaco_test{i},v0,lams,lame,ds0,'ds_max',ds_max,'display','off','l_0',(lame-lams)/2+lams);
    [probinfo,probcounter] = crashtest_check_output(vs,ls,exitflag,lams,lame,sprintf('### %d: l_0: (l_end-l_start)/2+l_start ###\n',i),probinfo,probcounter);
    %% l_target:
    fprintf('\n### %d: l_target: (l_end-l_start)/2+l_start ###\n',i);
    [vs,ls,exitflag] = continuation(fun_jaco_test{i},v0,lams,lame,ds0,'ds_max',ds_max,'display','off','l_target',(lame-lams)/2+lams);
    [probinfo,probcounter] = crashtest_check_output(vs,ls,exitflag,lams,lame,sprintf('### %d: l_target: (l_end-l_start)/2+l_start ###\n',i),probinfo,probcounter);
    %% l_target:
    fprintf('\n### %d: alpha_reverse: pi/4 ###\n',i);
    [vs,ls,exitflag] = continuation(fun_jaco_test{i},v0,lams,lame,ds0,'ds_max',ds_max,'display','off','alpha_reverse',pi/4);
    [probinfo,probcounter] = crashtest_check_output(vs,ls,exitflag,lams,lame,sprintf('### %d: alpha_reverse: pi/4 ###\n',i),probinfo,probcounter);
    %% plot:
    fprintf('\n### %d: plot: on ###\n',i);
    [vs,ls,exitflag] = continuation(fun_jaco_test{i},v0,lams,lame,ds0,'ds_max',ds_max,'display','off','plot','on');
    [probinfo,probcounter] = crashtest_check_output(vs,ls,exitflag,lams,lame,sprintf('### %d: plot: on ###\n',i),probinfo,probcounter);
    fprintf('\n### %d: plot: on & bifurcation: mark ###\n',i);
    [vs,ls,exitflag] = continuation(fun_jaco_test{i},v0,lams,lame,ds0,'ds_max',ds_max,'display','off','plot','on','bifurcation','mark');
    [probinfo,probcounter] = crashtest_check_output(vs,ls,exitflag,lams,lame,sprintf('### %d: plot: on & bifurcation: mark ###\n',i),probinfo,probcounter);
    fprintf('\n### %d: plot: on & bifurcation: determine ###\n',i);
    [vs,ls,exitflag] = continuation(fun_jaco_test{i},v0,lams,lame,ds0,'ds_max',ds_max,'display','off','plot','on','bifurcation','determine');
    [probinfo,probcounter] = crashtest_check_output(vs,ls,exitflag,lams,lame,sprintf('### %d: plot: on & bifurcation: determine ###\n',i),probinfo,probcounter);
    fprintf('\n### %d: plot: on & bifurcation: trace ###\n',i);
    [vs,ls,exitflag] = continuation(fun_jaco_test{i},v0,lams,lame,ds0,'ds_max',ds_max,'display','off','plot','on','bifurcation','trace');
    [probinfo,probcounter] = crashtest_check_output(vs,ls,exitflag,lams,lame,sprintf('### %d: plot: on & bifurcation: trace ###\n',i),probinfo,probcounter);
    fprintf('\n### %d: plot_vars_index: 1 ###\n',i);
    [vs,ls,exitflag] = continuation(fun_jaco_test{i},v0,lams,lame,ds0,'ds_max',ds_max,'display','off','plot','on','plot_vars_index',1);
    [probinfo,probcounter] = crashtest_check_output(vs,ls,exitflag,lams,lame,sprintf('### %d: plot_vars_index: 1 ###\n',i),probinfo,probcounter);
    %% include_reverse:
    fprintf('\n### %d: include_reverse: on ###\n',i);
    [vs,ls,exitflag] = continuation(fun_jaco_test{i},v0,lams,lame,ds0,'ds_max',ds_max,'display','off','include_reverse','on');
    [probinfo,probcounter] = crashtest_check_output(vs,ls,exitflag,lams,lame,sprintf('### %d: include_reverse: on ###\n',i),probinfo,probcounter);
    %% predictor_tangential:
    fprintf('\n### %d: predictor_tangential ###\n',i);
    [vs,ls,exitflag] = continuation(fun_jaco_test{i},v0,lams,lame,ds0,'ds_max',ds_max,'display','off','predictor','tangential');
    [probinfo,probcounter] = crashtest_check_output(vs,ls,exitflag,lams,lame,sprintf('### %d: predictor_tangential ###\n',i),probinfo,probcounter);
    %% predictor_polynomial_order:
    fprintf('\n### %d: predictor_polynomial_order: 2 ###\n',i);
    [vs,ls,exitflag] = continuation(fun_jaco_test{i},v0,lams,lame,ds0,'ds_max',ds_max,'display','off','predictor','polynomial','predictor_polynomial_order',2);
    [probinfo,probcounter] = crashtest_check_output(vs,ls,exitflag,lams,lame,sprintf('### %d: predictor_polynomial_order: 2 ###\n',i),probinfo,probcounter);
    %% predictor_polynomial_fit:
    fprintf('\n### %d: predictor_polynomial_fit: 4 ###\n',i);
    [vs,ls,exitflag] = continuation(fun_jaco_test{i},v0,lams,lame,ds0,'ds_max',ds_max,'display','off','predictor','polynomial','predictor_polynomial_fit',4);
    [probinfo,probcounter] = crashtest_check_output(vs,ls,exitflag,lams,lame,sprintf('### %d: predictor_polynomial_fit: 4 ###\n',i),probinfo,probcounter);
    %% predictor_polynomial_adaptive:
    fprintf('\n### %d: predictor_polynomial_adaptive: nt = 4, nf = 5 ###\n',i);
    [vs,ls,exitflag] = continuation(fun_jaco_test{i},v0,lams,lame,ds0,'ds_max',ds_max,'display','off','predictor','polynomial','predictor_polynomial_adaptive','on','predictor_polynomial_order',4,'predictor_polynomial_fit',5);
    [probinfo,probcounter] = crashtest_check_output(vs,ls,exitflag,lams,lame,sprintf('### %d: predictor_polynomial_adaptive: nt = 4, nf = 5 ###\n',i),probinfo,probcounter);
    %% predictor_polynomial_adaptive:
    fprintf('\n### %d: predictor_solver: on ###\n',i);
    [vs,ls,exitflag] = continuation(fun_jaco_test{i},v0,lams,lame,ds0,'ds_max',ds_max,'display','off','predictor_solver','on');
    [probinfo,probcounter] = crashtest_check_output(vs,ls,exitflag,lams,lame,sprintf('### %d: predictor_solver: on ###\n',i),probinfo,probcounter);
    %% step_size_control:
    fprintf('\n### %d: step_size_control: angle ###\n',i);
    [vs,ls,exitflag] = continuation(fun_jaco_test{i},v0,lams,lame,ds0,'ds_max',ds_max,'display','off','step_size_control','angle');
    [probinfo,probcounter] = crashtest_check_output(vs,ls,exitflag,lams,lame,sprintf('### %d: step_size_control: angle ###\n',i),probinfo,probcounter);
    fprintf('\n### %d: step_size_control: pid ###\n',i);
    [vs,ls,exitflag] = continuation(fun_jaco_test{i},v0,lams,lame,ds0,'ds_max',ds_max,'display','off','step_size_control','pid');
    [probinfo,probcounter] = crashtest_check_output(vs,ls,exitflag,lams,lame,sprintf('### %d: step_size_control: pid ###\n',i),probinfo,probcounter);
    fprintf('\n### %d: step_size_control: curvature ###\n',i);
    [vs,ls,exitflag] = continuation(fun_jaco_test{i},v0,lams,lame,ds0,'ds_max',ds_max,'display','off','step_size_control','curvature');
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