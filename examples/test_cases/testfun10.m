addpath('test_cases\testfun10_aux');

v0 = 10^-15;
lams = 1;
lame = 10;
ds0 = 0.05;
ds_max = 0.1;

opt = optimoptions('fsolve','display','None');
solver = @(f,x) fsolve(f,x,opt);

fun = @(v,l) residual_fun10(v,l,solver);