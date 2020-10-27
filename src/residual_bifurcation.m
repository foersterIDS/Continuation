%% path continuation - residual_bifurcation
%
%   Institute of Dynamics and Vibration Research
%   Leibniz University Hannover
%   27.10.2020 - Alwin Förster
%
function [R] = residual_bifurcation(fun,x,Opt,scale)
    if Opt.jacobian
        [R1,J1] = fun(x(1:end-1),x(end));
        if diff(size(J1))
            J1 = J1(:,1:end-1);
        end
    else
        R1 = fun(x(1:end-1),x(end));
        fun_J = @(v) fun(v,x(end));
        J1 = numeric_jacobian(fun_J,x(1:end-1),R1);
    end
    R2 = det(J1)*scale;
    R = [R1;R2];
end