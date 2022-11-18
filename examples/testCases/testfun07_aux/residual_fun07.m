function [f,J] = residual_fun07(v,l,n)
%     5-2*(x-5)^2+(x-5)^3
    a = 5;
    b = -2;
    c = 1;
    d = -5;
    sc = ones(n,1);%(2+sin(1:n))';
    f = sc*l-[a+b*(v+d).^2+c*(v+d).^3];
    J = [diag(-((d+v).*(2*b+3*c*(d+v)))),sc];
end