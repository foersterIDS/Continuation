addpath('test_cases\testfun03_aux');
%addpath('..\..');
%addpath('..\..\ELM');
%addpath('..\..\ANSYS');
%addpath('..\..\pathConti');
%addpath('..\..\Fixpunkt');

is = [];
ts = [];
iter_max = 10^4;
duf = 'DUF';
ecs = 'ECS';
nle_typ = ecs; % 'DUF' or 'ECS'

%% Settings:
omu = 0.5;
omo = 2;%1.6;
lams = omu;
lame = omo;

%% System:
Ndof = 1;
nex = Ndof;
if strcmp(nle_typ,duf)
    nx = 2*Ndof;
elseif strcmp(nle_typ,ecs)
    nx = 2*Ndof+1;
end
nz = nx+2;

m = 1;%1.2;
c = 0.05;%0.9;
k = 1;%580.8;
kt = 2;%230.4;
l = 0.1;%kt;
kp = 0;
mu = 0.3;
fN = 2;%2;%100;

Df = 0.01;
sf = 0.18/sqrt(2);%15.23;
ab = [0,1];
om = omu;%sqrt((k+kp)/m);

M = m*eye(Ndof);
K = kette(k+kp,Ndof,1);
C = c/(k+kp)*K;

fi = Filter2(nex,Df,om,sf,ab);
fex = FilteredProcess(zeros(Ndof,1),fi);
SM2O = eye(Ndof);

[Ax,Bx,SM1Ox,dim_infox] = NonLinearSystem_1O.getAB(M,C,K,SM2O);

fnlx = NonLinearForce(nx);
if strcmp(nle_typ,duf)
    nle = NL_Element_DUF(l);
    fnlx.addElement(nle,Ndof+1,1);
    nls = NonLinearSystem_1O(Ax,Bx,fnlx,SM1Ox,fex,dim_infox);
    nls.setMeanFree(false);
elseif strcmp(nle_typ,ecs)
    nle = NL_Element_ECS(mu,fN,kt);
    fnlx.addElement(nle,[Ndof+1,nx],[Ndof+1,nx]);
    nls = NonLinearSystem_1O(Ax,Bx,fnlx,SM1Ox,fex,dim_infox);
    nls.setMeanFree(true);
end

f = @(vars) fixedpointfun_ELM(nls,vars);
fc = @(vars) constrainFun_ELM(nls,vars);
x0 = covToVars(nls,eye(nls.nx),0*ones(nls.nx,1));%[5;0.03;0.01;0;0];%

try
    x = fixedpoint(f,x0,'constraint',fc,'damping','off','tol',10^-8,'itermax',20);
catch
    x = x0;
end

[Kzz_red,muz_red] = varsToCov(nls,x);
ls = nls.getLinearizedSystem(muz_red,Kzz_red);
Kzz = ls.getKXX();

Kzz(nls.dim_info==0,nls.dim_info==0)

%% pathConti
fun = @(vars,om) residual_filter2(nls.nlsh,fex,vars,om);
ds0 = 0.0051;
v0 = x;