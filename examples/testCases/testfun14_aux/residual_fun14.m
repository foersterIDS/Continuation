function [varargout] = residual_fun14(v,l,a,is_jacobian)
    Res = [v(1) - a*sin(2*pi*l);
           v(2) - a*sin(2*pi*l)*cos(2*pi*l)];
    varargout{1} = Res;
    if is_jacobian
        J = [1, 0, -2*pi*a*cos(2*pi*l);
             0, 1, -2*pi*a*cos(4*pi*l)];
        varargout{2} = J;
    end
end