function [varargout] = residual_fun05(v,l,is_jacobian,r,lm,vm)
    f = ((l-lm)^2+(v-vm)^2-r^2)*(v-l^2);
    varargout{1} = f;
    if is_jacobian
        J = [(2*(v-vm))*(v-l^2)+((l-lm)^2+(v-vm)^2-r^2)*(1),(2*(l-lm))*(v-l^2)-((l-lm)^2+(v-vm)^2-r^2)*(2*l)];
        varargout{2} = J;
    end
end