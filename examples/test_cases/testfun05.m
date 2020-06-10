addpath('test_cases\testfun05_aux');
lm = 2;
vm = lm^3;
r = lm/2;
fun = @(v,l) residual_fun05(v,l,true,r,lm,vm);
lams = 0;
lame = 2*lm;
v0 = 0;
ds0 = 0.05;
ds_max = 0.1;
%% Test-function with and without jacobian
fun_jaco_test = cell(2,1);
fun_jaco_test{1} = @(v,l) residual_fun05(v,l,true,r,lm,vm);
fun_jaco_test{2} = @(v,l) residual_fun05(v,l,false,r,lm,vm);