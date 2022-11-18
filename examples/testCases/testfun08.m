addpath('testCases\testfun08_aux');
r = 1;
fun = @(v,l) residual_fun08(v,l,true,r);
lams = 0;
lame = 2*r;
v0 = r/2;
ds0 = r/100;
dsMax = r/10;
%% Test-function with and without jacobian
funJacoTest = cell(2,1);
funJacoTest{1} = @(v,l) residual_fun08(v,l,true,r);
funJacoTest{2} = @(v,l) residual_fun08(v,l,false,r);