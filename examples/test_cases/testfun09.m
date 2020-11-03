addpath('test_cases\testfun09_aux');
vm = 0;
lm = 0;
r = 1;
v0 = 1.02;
lams = -1;
lame = 4;
ds0 = 0.05;
ds_max = 0.1;

fun = @(v,l) residual_fun09(v,l,r,vm,lm);