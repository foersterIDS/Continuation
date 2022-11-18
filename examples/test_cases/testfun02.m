%% one-dimensional duffing oscillator
% mu \ddot q + zeta \dot q + kappa q + \gamma q^3 = P cos( l * t )
%
% 26.07.2022      Anna Lefken

addpath('test_cases\testfun02_aux');
%% Parameters
mu = 1;
zeta = 0.1;
kappa = 1;
gamma = 0.5;
P=4;

Nh=5;                                                   % number of harmonics
Np=2^8;                                                 % number of sampling points per period
lams = 0;                                             % lambda at start of continuation
lame = 6;                                             % lambda at end of continuation
ds0 = mu/100;                                           % start step length
dsMax = mu/10;                                         % maximum step length
%% Calculation
[G, H] = func_FourierMatrix(Nh,Np);

p=zeros(2*(Nh)+1,1);                                    % p vector in frequency domain
p(2)=P;

A1=lams*kron(diag(1:Nh),[0 1;-1 0]);
A=[0 zeros(1, 2*Nh);zeros(2*Nh,1) A1];
Sh=kron(A^2,mu)+kron(A,zeta)+kron(A^0,kappa);           % dynamic stiffness at start

v0=Sh\p;                                                % start values from linear system
%% Function
fun = @(v,l) residual_fun02(v,l,mu,zeta,kappa,gamma,p,Nh,G,H);
