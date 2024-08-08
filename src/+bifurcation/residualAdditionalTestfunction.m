%% path continuation - bifurcation.residualAdditionalTestfunction
%
%   Institute of Dynamics and Vibration Research
%   Leibniz University Hannover
%   25.05.2023 - Anna Lefken
%
function [R] = residualAdditionalTestfunction(func,x,Opt,Path,Info)
    if Opt.jacobian
        [R1,~] = func(x(1:Info.nv),x(Info.nv+1));
    else
        R1 = func(x(1:end-1),x(end));
        
    end
     % evaluate Jacobian for every new x
     % cannot use input "Jacobian", as this is evaluated at the last path step
    Jacobianbif.solver = aux.numericJacobian(@(x)func(x(1:end-1),x(end)), x,'diffquot', Opt.diffquot);                  
    R2 = Opt.bifAdditionalTestfunction(func,x,Jacobianbif,Path,Info);
    R = [R1;R2];
end