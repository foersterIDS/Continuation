addpath('test_cases\testfun04_aux');
fun = @(v,l) residual_fun04(v,l,true);
lams = 1;
lame = 10;
v0 = [0;1];
ds0 = 1;
dsMax = 2;
%% Test-function with and without jacobian
fun_jaco_test = cell(2,1);
fun_jaco_test{1} = @(v,l) residual_fun04(v,l,true);
fun_jaco_test{2} = @(v,l) residual_fun04(v,l,false);