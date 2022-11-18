addpath('test_cases\testfun14_aux');

a = 5;
is_jacobian = 1;
lams = 0;
lame = 10;
v0 = [0,0].';
ds0 = 0.3;
dsMax = 0.5;


fun = @(v,l) residual_fun14(v,l,a,is_jacobian);


