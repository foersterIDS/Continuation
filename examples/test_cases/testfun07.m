addpath('test_cases\testfun07_aux');
n = 500;
fun = @(v,l) residual_fun07(v,l,n);
lams = 0;
lame = 10;
v0 = fsolve(@(v) fun(v,lams),3.75*ones(n,1),optimoptions('fsolve','display','off'));
ds0 = 0.01;
dsMax = 1;