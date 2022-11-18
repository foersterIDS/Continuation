%% Residuum for one-dimensional duffing oscillator
function R=residual_fun02(v,l,mu,zeta,kappa,gamma,p,Nh,G,H)

A1=l*kron(diag(1:Nh),[0 1;-1 0]);
A=[0 zeros(1, 2*Nh);zeros(2*Nh,1) A1];
Sh=kron(A^2,mu)+kron(A,zeta)+kron(A^0,kappa);       % dynamic stiffness matrix

v_tau=v'*H;                                         % transform to time domain
fnl_tau=gamma*v_tau.^3;                             % force of cubic stiffness
fnl=G*fnl_tau';                                     % transform to frequency domain

%% Residuum
R=Sh*v+fnl-p;


end