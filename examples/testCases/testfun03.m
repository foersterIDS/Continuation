addpath('testCases/testfun03_aux');

%% Settings:
omu = 0.5;
omo = 1.6;
lams = omu;
lame = omo;

%% System:
Ndof = 1;
nex = Ndof;
nx = 2*Ndof;
nz = nx+2;

m = 1;
c = 0.05;
k = 1;
l = 0.1;

Df = 0.01;
sf = 0.18/sqrt(2);
om = omu;

Az = [1,0,0,0;
      0,m,0,0;
      0,0,1,0;
      0,0,0,1];
invAz = [  1,  0,  0,  0;
           0,1/m,  0,  0;
           0,  0,  1,  0;
           0,  0,  0,  1];
Bz = @(v,w) [      0,     -1,      0,      0;
             k+3*l*v,      c,      0,      1;
                   0,      0,      0,     -1;
                   0,      0,    w^2, 2*Df*w];
Bz_dpa = @(v,w,l) [      0,     -1,      0,      0;
                   k+3*l*v,      c,      0,      1;
                         0,      0,      0,     -1;
                         0,      0,    w^2, 2*Df*w];
SMz = [0;0;0;1];
Sff = @(w) 4*Df*w*sf^2;

G = @(v,w) -invAz*Bz(v,w);
G_dpa = @(v,w,l) -invAz*Bz_dpa(v,w,l);
D = @(v,w) conj(invAz*SMz)*Sff(w)*(invAz*SMz).';

vs = [1;0;0;0];

%% pathConti
fun = @(vars,om) vs'*lyap(G(vars,om),D(vars,om))*vs-vars;
fun_dpa = @(vars,om,g) vs'*lyap(G_dpa(vars,om,g),D(vars,om))*vs-vars;
ds0 = 0.0001;
dsMax = 0.01;
v0 = 0.05;