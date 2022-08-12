%% path continuation - bifurcation.residual
%
%   Institute of Dynamics and Vibration Research
%   Leibniz University Hannover
%   27.10.2020 - Alwin Förster
%
function [R] = residual_additional_testfunction(func,x,Opt,Jacobian,Path,Info)
    if Opt.jacobian
        [R1,~] = func(x(1:Info.nv),x(Info.nv+1));
    else
        R1 = func(x(1:end-1),x(end));
        
    end
    R2 = Opt.bif_additional_testfunction(func,x,Jacobian,Path,Info);
    R = [R1;R2];
end