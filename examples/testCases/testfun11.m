addpath('testCases/testfun11_aux');

v0 = 10^-15;
lams = 0;
lame = 10;
ds0 = 0.05;
dsMax = 0.1;

fun = @(v,l) residual_fun11(v,l);