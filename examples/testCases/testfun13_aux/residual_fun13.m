function [R,J] = residual_fun13(v,l)
    R = v-sin(1/l);
    J = [1,cos(1/l)/l^2];
end