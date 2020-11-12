function [varargout] = residual_fun08(v,l,is_jacobian,r)
    f = (v-r/2-1/2*l).*(2*r-2*v+1/r*v.^2-l);
    varargout{1} = f;
    if is_jacobian
        J = [(2*r-2*v+1/r*v.^2-l)+(v-r/2-1/2*l).*(-2+2*1/r*v),(-1/2).*(2*r-2*v+1/r*v.^2-l)+(v-r/2-1/2*l).*(-1)];
        varargout{2} = J;
    end
end