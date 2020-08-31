addpath('test_cases\testfun06_aux');
Nres = 5;
fun = @(v,l) residual_fun06(v,l,Nres);
lams = 0.1;
lame = 1;
v0 = ones(Nres,1);
ds0 = 0.01;
ds_max = 0.1;