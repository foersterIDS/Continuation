%% path continuation - bifurcation.residual
%
%   Institute of Dynamics and Vibration Research
%   Leibniz University Hannover
%   27.10.2020 - Alwin Förster
%
function [R] = residual(fun,x,oih)
    if oih.opt.jacobian
        [R1,J1] = fun(x(1:oih.info.nv),x(oih.info.nv+1));
        if diff(size(J1))
            J1 = J1(1:oih.info.nv,1:oih.info.nv);
        end
    else
        R1 = fun(x(1:end-1),x(end));
        funJ = @(v) fun(v,x(end));
        J1 = aux.numericJacobian(funJ,x(1:oih.info.nv),'diffquot',oih.opt.diffquot,'diffStep',oih.opt.diffStep);
    end
    
    if oih.opt.bifResidual.determinant
        R2 = det(J1)*oih.bifurcation.scaling(end);
    elseif oih.opt.bifResidual.luFactorization
        [~,U] = lu(J1);
        diagU=diag(U);
        [~,indMin]=min(abs(diagU));
        R2 = diagU(indMin);
    end
    R = [R1;R2];
end