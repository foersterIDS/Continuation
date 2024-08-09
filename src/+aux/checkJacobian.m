%% path continuation - aux.checkJacobian
%
%   Institute of Dynamics and Vibration Research
%   Leibniz University Hannover
%   25.01.2024 - Alwin Förster
%
function checkJacobian(fun,f0,x0,oih)
    %% arguments
    arguments
        fun (1,1) function_handle
        f0 (:,1) double
        x0 (:,1) double
        oih (1,1) aux.OptInfoHandle
    end
    %% check
    if oih.opt.jacobian && oih.opt.checkJacobian
        %% check validity
        if isMATLABReleaseOlderThan('R2023b')
            %% R2023a and older
            % numeric Jacobian
            jacobianNum = aux.numericJacobian(@(xt) fun(xt),x0,'centralValue',f0,'diffquot',oih.opt.diffquot);
            % compare with user provided Jacobian
            [na1,na2] = size(oih.solver.jacobian);
            [nn1,nn2] = size(jacobianNum);
            n1 = min(na1,nn1);
            n2 = min(na2,nn2);
            dJrel = abs(oih.solver.jacobian(1:n1,1:n2)-jacobianNum(1:n1,1:n2))./round(abs(oih.solver.jacobian(1:n1,1:n2)),8);
            err = norm(dJrel(logical((~isinf(dJrel)).*(~isnan(dJrel)))));
            valid = err<10^-3;
        else
            %% since R2023b
            [valid,temp] = checkGradients(fun,x0,'Tolerance',10^-4);
            err = temp.Objective;
        end
        %% set info & solver
        if valid && ~oih.info.validJacobian
            aux.printLine(oih,'--> Valid user provided Jacobian. Using user provided Jacobian for further computations.\n');
        elseif ~valid && oih.info.validJacobian
            aux.printLine(oih,'--> Invalid user provided Jacobian. Using numeric Jacobian for further computations.\n');
        end
        if valid
            oih.info.validJacobian = true;
            if ~isempty(oih.solver.temp)
                oih.solver.main = oih.solver.temp;
                oih.solver.temp = [];
            end
        else
            oih.info.validJacobian = false;
            oih.solver.temp = oih.solver.main;
            oih.solver.main = oih.solver.numJac;
        end
    end
end