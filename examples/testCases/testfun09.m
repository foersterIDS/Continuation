addpath('testCases\testfun09_aux');
vm = 0;
lm = 0;
r = 1;
v0 = 0.001;
lams = -1;
lame = 4;
ds0 = 0.05;
dsMax = 0.1;

fun = @(v,l) residual_fun09(v,l,r,vm,lm);
% benoetigt 'direction',[1;0]