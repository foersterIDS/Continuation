function [R,J] = res_fun(xsc,sc)
    x = xsc.*sc;
    r1 = 10^+3;
    r2 = 10^-3;
    x2 = r2/2;
    R = [(x(1)/r1)^2+(x(2)/r2)^2-1;x(2)-x2];
    J = [2*x(1)/r1^2,2*x(2)/r2^2;0,1];
    J = J*diag(sc);
end