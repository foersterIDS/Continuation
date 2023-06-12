addpath('testCases/testfun12_aux');

v0 = [-1;0;0];
% v0 = -1;
lams = -3;
lame = 10;
ds0 = 0.05;
dsMax = 0.1;

fun = @(v,l) residual_fun12(v,l);