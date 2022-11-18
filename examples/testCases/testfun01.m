addpath('testCases\testfun01_aux');
fun = @(v,l) residual_fun01(v,l);
lams = 0.1;
lame = 1;
v0 = 4;
ds0 = 0.01;
dsMax = 0.1;