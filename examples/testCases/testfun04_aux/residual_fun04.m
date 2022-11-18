function [varargout] = residual_fun04(v,l,is_jacobian)
    f = [v'*v-l^2;v(2)-sin(l)*exp(v(1))];
    varargout{1} = f;
    if is_jacobian
        J = [2*v',-2*l;-sin(l)*exp(v(1)),1,-cos(l)*exp(v(1))];
        varargout{2} = J;
    end
end