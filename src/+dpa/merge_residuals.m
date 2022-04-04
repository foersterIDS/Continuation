%% path continuation - dpa.merge_residuals
%
%   Institute of Dynamics and Vibration Research
%   Leibniz University Hannover
%   31.03.2022 - Alwin Förster
%
function [varargout] = merge_residuals(Opt,fun,res_dpa,x,g)
    nv = numel(x)-1;
    v = x(1:nv);
    l = x(nv+1);
    z = [x;g];
    if Opt.jacobian
        [Rv,Jfi] = fun(v,l,g);
        [nc,ny] = size(Jfi);
        if nc~=nv
            error('Jacobian has invalid number of rows.');
        end
        funJ = @(z) fun(z(1:nv),z(nv+1),z(nv+2));
        Jfi_add = aux.numeric_jacobian(funJ,z,'derivative_dimensions',(ny+1):(nv+2));
        Jfz = [Jfi,Jfi_add];
        funl = @(z) res_dpa(z(1:nv),z(nv+1),z(nv+2));
        Rl = funl(z);
        Jlz = aux.numeric_jacobian(funl,z,'derivative_dimensions',1:(nv+2));
        varargout{1} = [Rv;Rl];
        varargout{2} = [Jfz;Jlz];
    else
        Rv = fun(v,l,g);
        Rl = res_dpa(v,l,g);
        varargout{1} = [Rv;Rl];
    end
end