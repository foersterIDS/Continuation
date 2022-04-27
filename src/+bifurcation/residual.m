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
        fun_J = @(v) fun(v,x(end));
        J1 = aux.numeric_jacobian(fun_J,x(1:Info.nv),R1,'diffquot',Opt.diffquot);
    end
    R2 = det(J1)*scale;
    R = [R1;R2];
end