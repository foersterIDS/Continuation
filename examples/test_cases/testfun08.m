addpath('test_cases\testfun08_aux');
r = 1;
fun = @(v,l) residual_fun08(v,l,r);
lams = 0;
lame = 2*r;
v0 = r/2;
ds0 = r/100;
ds_max = r/10;