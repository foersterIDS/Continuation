%% path continuation - dpa.mergeResiduals
%
%   Institute of Dynamics and Vibration Research
%   Leibniz University Hannover
%   31.03.2022 - Alwin Förster
%
function [varargout] = mergeResiduals(oih,fun,resDpa,x,g)
    nv = numel(x)-1;
    v = x(1:nv);
    l = x(nv+1);
    z = [x;g];
    if oih.opt.jacobian
        [Rv,Jfi] = fun(v,l,g);
        [nc,ny] = size(Jfi);
        if nc~=nv
            error('Jacobian has invalid number of rows.');
        end
        funJ = @(z) fun(z(1:nv),z(nv+1),z(nv+2));
        JfiAdd = aux.numericJacobian(funJ,z,'derivativeDimensions',(ny+1):(nv+2),'diffStep',oih.opt.diffStep);
        Jfz = [Jfi,JfiAdd];
        funl = @(z) resDpa(z(1:nv),z(nv+1),z(nv+2));
        Rl = funl(z);
        Jlz = aux.numericJacobian(funl,z,'derivativeDimensions',1:(nv+2),'diffStep',oih.opt.diffStep);
        varargout{1} = [Rv;Rl];
        varargout{2} = [Jfz;Jlz];
    else
        Rv = fun(v,l,g);
        Rl = resDpa(v,l,g);
        varargout{1} = [Rv;Rl];
    end
end