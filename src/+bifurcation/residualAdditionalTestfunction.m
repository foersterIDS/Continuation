%% path continuation - bifurcation.residualAdditionalTestfunction
%
%   Institute of Dynamics and Vibration Research
%   Leibniz University Hannover
%   25.05.2023 - Anna Lefken
%
function [R] = residualAdditionalTestfunction(func,x,oih)
    if oih.opt.jacobian
        [R1,~] = func(x(1:oih.info.nv),x(oih.info.nv+1));
    else
        R1 = func(x(1:end-1),x(end));
        
    end
     % evaluate Jacobian for every new x
     % cannot use input "Jacobian", as this is evaluated at the last path step
    Jacobianbif.solver = aux.numericJacobian(@(x)func(x(1:end-1),x(end)), x,'diffquot', oih.opt.diffquot);                  
    R2 = oih.opt.bifAdditionalTestfunction(func,x,Jacobianbif,oih.path,oih.info);
    R = [R1;R2];
end