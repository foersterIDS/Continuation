%% path continuation - bifurcation.residual
%
%   Institute of Dynamics and Vibration Research
%   Leibniz University Hannover
%   27.10.2020 - Alwin Förster
%
function [R] = residual(fun,x,Opt,Info,scale)
    if Opt.jacobian
        [R1,J1] = fun(x(1:Info.nv),x(Info.nv+1));
        if diff(size(J1))
            J1 = J1(1:Info.nv,1:Info.nv);
        end
    else
        R1 = fun(x(1:end-1),x(end));
        funJ = @(v) fun(v,x(end));
        J1 = aux.numericJacobian(funJ,x(1:Info.nv),R1,'diffquot',Opt.diffquot);
    end
    
    if Opt.bifResidual.determinant
        R2 = det(J1)*scale;
    elseif Opt.bifResidual.luFactorization
        [~,U] = lu(J1);
        diagU=diag(U);
        [~,indMin]=min(abs(diagU));
        R2 = diagU(indMin);
    end
    R = [R1;R2];
end