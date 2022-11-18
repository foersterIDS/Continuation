function [f,J] = residual_fun01(v,l)
    f = v.^2+5-exp((1/l)*v);
    J = [2*v-1/l*exp(v/l),v/l^2*exp(v/l)];
end