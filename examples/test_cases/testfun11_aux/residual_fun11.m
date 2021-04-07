function [R,J] = residual_fun11(v,l)
    R = v-10-10^-10*sin(2*pi/10*l);
    J = [1,-10^-10*2*pi/10*cos(2*pi/10*l)];
end