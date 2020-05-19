addpath('test_cases\testfun02_aux');
mu = 1;%0.5;%1;
zeta = 0.05;%0.3;%0.05;
kappa = 1;%0.7;%1;
gamma = 0.1;%5;%0.1;
P = 0.18;%1;%0.18;
H = 7;      % harmonic order
N = 2^6;    % number of time samples per period
lams = 0.5;%0.001;%0.5;  % start frequency
lame = 1.6;%4;%1.6; % end frequency
Q = (-lams^2*mu+1i*lams*zeta+kappa)\P;
v0 = [0;real(Q);-imag(Q);zeros(2*(H-1),1)];
ds0 = 0.1;
fun = @(x,l) HB_residual_Duffing([x;l],mu,zeta,kappa,gamma,P,H,N);