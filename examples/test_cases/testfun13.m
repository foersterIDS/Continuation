addpath('test_cases\testfun13_aux');

lams = 100;
lame = 0;
v0 = sin(1/lams);
ds0 = 0.05;
dsMax = 0.1;

fun = @(v,l) residual_fun13(v,l);