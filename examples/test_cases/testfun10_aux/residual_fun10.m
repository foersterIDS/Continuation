function [R,J] = residual_fun10(v,l,solver)
    %% set parameters:
    damper = 0.0;
    m = 1;
    k = 15;
    c = 0.1;
    knl = 78;
    vz1 = +1;
    vz2 = -1;
    du1 = +1;%*10^-3;
    du2 = -du1;
    sf = 1;%10^-3;%1
    omf = l;
    Df = 0.01;
    SMq = 1;
    a = 0;
    b = 1;
    S0 = (4*Df*omf*sf^2)/(a^2/omf^2+b^2);
    dS0dl = (4*Df*sf^2*omf^2*(3*a^2+b^2*omf^2))/(a^2+b^2*omf^2)^2;
    Kqq = v;
    %% expand:
    SMx = [0;SMq];
    SMy = [0;1];
    SMz = [0;0;SMy];
    Ax = m*eye(2);
    Bx = [0,-m;k,c];
    Ay = eye(2);
    By = [0,-1;omf^2,2*Df*omf];
    Az = [Ax,zeros(2);zeros(2),Ay];
    invAz = inv(Az);
    Bz = [Bx,-a*SMx,-b*SMx;zeros(2),By];
    Bnlz = zeros(4);
    mufnlq = @(muq) (knl*(vz1*sqrt(Kqq/(2*pi))*exp(-(muq-du1)^2/(2*Kqq))+(muq-du1)/2*(1+vz1*erf((muq-du1)/sqrt(2*Kqq)))))+...
                    (knl*(vz2*sqrt(Kqq/(2*pi))*exp(-(muq-du2)^2/(2*Kqq))+(muq-du2)/2*(1+vz2*erf((muq-du2)/sqrt(2*Kqq)))));
    fun_muq = @(muq) k*muq+mufnlq(muq);
    muq = solver(fun_muq,0);
    Bnlz(2,1) = (knl/2*(1+vz1*erf((muq-du1)/sqrt(2*Kqq)))+knl*vz1*du1/sqrt(2*pi*Kqq)*exp(-(muq-du1)^2/(2*Kqq)))+...
                (knl/2*(1+vz2*erf((muq-du2)/sqrt(2*Kqq)))+knl*vz2*du2/sqrt(2*pi*Kqq)*exp(-(muq-du2)^2/(2*Kqq)));
    %% lyap:
    Gz = -invAz*(Bz+Bnlz);
    Dz = conj(invAz*SMz)*S0*(invAz*SMz).';
    Kzz = lyap(Gz,Dz);
    Kqqip1 = Kzz(1,1);
    ffix = damper*Kqq+(1-damper)*Kqqip1;
    R = ffix-Kqq;
    %% Jacobi:
    dmufnlqdv = @(dmuqdv) (knl*vz1/sqrt(2*pi*Kqq)*exp(-(muq-du1)^2/(2*Kqq))*(1/2-(1/(2*Kqq))*(muq-du1)*(2*Kqq*dmuqdv-muq+du1)+(muq-du1)*(dmuqdv-(muq-du1)/(2*Kqq)))+knl/2*dmuqdv*(1+vz1*erf((muq-du1)/sqrt(2*Kqq))))+...
                          (knl*vz2/sqrt(2*pi*Kqq)*exp(-(muq-du2)^2/(2*Kqq))*(1/2-(1/(2*Kqq))*(muq-du2)*(2*Kqq*dmuqdv-muq+du2)+(muq-du2)*(dmuqdv-(muq-du2)/(2*Kqq)))+knl/2*dmuqdv*(1+vz2*erf((muq-du2)/sqrt(2*Kqq))));
    fun_dmuqdv = @(dmuqdv) k*dmuqdv+dmufnlqdv(dmuqdv);
    dmuqdv = solver(fun_dmuqdv,0);
    dBnlz21dv = (knl*vz1/sqrt(2*pi*Kqq)*exp(-(muq-du1)^2/(2*Kqq))*(dmuqdv-(muq-du1)/(2*Kqq)+du1*(-1/(2*Kqq)+(muq-du1)^2/(2*Kqq^2)-(muq-du1)*dmuqdv/Kqq)))+...
                (knl*vz2/sqrt(2*pi*Kqq)*exp(-(muq-du2)^2/(2*Kqq))*(dmuqdv-(muq-du2)/(2*Kqq)+du2*(-1/(2*Kqq)+(muq-du2)^2/(2*Kqq^2)-(muq-du2)*dmuqdv/Kqq)));
    dBnlzdv = zeros(4);
    dBnlzdv(2,1) = dBnlz21dv;
    dBnlzdl = zeros(4);
    dBnlzdl(4,3) = 2*omf;
    dBnlzdl(4,4) = 2*Df;
    dGzdv = -invAz*dBnlzdv;
    dGzdl = -invAz*dBnlzdl;
    dKzzdv = lyap(Gz,dGzdv*Kzz+Kzz*dGzdv');
    dKzzdl = lyap(Gz,dGzdl*Kzz+Kzz*dGzdl');
    dffixdv = damper+(1-damper)*dKzzdv(1,1);
    dffixdl = (1-damper)*dKzzdl(1,1);
    J = [dffixdv-1,dffixdl];
end