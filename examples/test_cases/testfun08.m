addpath('test_cases\testfun08_aux');
r = 1;
fun = @(v,l) residual_fun08(v,l,true,r);
lams = 0;
lame = 2*r;
v0 = r/2;
ds0 = r/100;
ds_max = r/10;
%% Test-function with and without jacobian
fun_jaco_test = cell(2,1);
fun_jaco_test{1} = @(v,l) residual_fun08(v,l,true,r);
fun_jaco_test{2} = @(v,l) residual_fun08(v,l,false,r);