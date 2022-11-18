function [f,J] = residual_fun09(v,l,r,vm,lm)
    f = (v-vm)^2 + (l-lm)^2 - r;
    J = [2*(v-vm), 2*(l-lm)];
end