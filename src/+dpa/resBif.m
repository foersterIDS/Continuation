%% path continuation - dpa.resBif
%
%   Institute of Dynamics and Vibration Research
%   Leibniz University Hannover
%   31.03.2022 - Alwin Förster
%
function [R] = resBif(fun,x,g,Opt,scale)
    if Opt.jacobian
        [~,J1] = fun(x(1:end-1),x(end),g);
        if diff(size(J1))
            J1 = J1(:,1:end-1);
        end
    else
        R1 = fun(x(1:end-1),x(end),g);
        funJ = @(v) fun(v,x(end),g);
        J1 = aux.numericJacobian(funJ,x(1:end-1),R1,'diffquot',Opt.diffquot);
    end
    R = det(J1)*scale;
end