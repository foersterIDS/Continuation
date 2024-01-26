%% path continuation - aux.checkJacobian
%
%   Institute of Dynamics and Vibration Research
%   Leibniz University Hannover
%   25.01.2024 - Alwin Förster
%
function [Info,Solver] = checkJacobian(fun,f0,x0,Info,Jacobian,Opt,Solver)
    %% arguments
    arguments
        fun (1,1) function_handle
        f0 (:,1) double
        x0 (:,1) double
        Info (1,1) struct
        Jacobian (1,1) struct
        Opt (1,1) struct
        Solver (1,1) struct
    end
    %% check
    if Opt.jacobian && Opt.checkJacobian
        %% check validity
        if isMATLABReleaseOlderThan('R2023b')
            %% R2023a and older
            % numeric Jacobian
            jacobianNum = aux.numericJacobian(@(xt) fun(xt),x0,'centralValue',f0,'diffquot',Opt.diffquot);
            % compare with user provided Jacobian
            [na1,na2] = size(Jacobian.solver);
            [nn1,nn2] = size(jacobianNum);
            n1 = min(na1,nn1);
            n2 = min(na2,nn2);
            dJrel = abs(Jacobian.solver(1:n1,1:n2)-jacobianNum(1:n1,1:n2))./round(abs(Jacobian.solver(1:n1,1:n2)),8);
            err = norm(dJrel(logical((~isinf(dJrel)).*(~isnan(dJrel)))));
            valid = err<10^-3;
        else
            %% since R2023b
            [valid,temp] = checkGradients(fun,x0,'Tolerance',10^-4);
            err = temp.Objective;
        end
        %% set info & solver
        if valid && ~Info.validJacobian
            aux.printLine(Opt,'--> Valid user provided Jacobian. Using user provided Jacobian for further computations.\n');
        elseif ~valid && Info.validJacobian
            aux.printLine(Opt,'--> Invalid user provided Jacobian. Using numeric Jacobian for further computations.\n');
        end
        if valid
            Info.validJacobian = true;
            if ~isempty(Solver.temp)
                Solver.main = Solver.temp;
                Solver.temp = [];
            end
        else
            Info.validJacobian = false;
            Solver.temp = Solver.main;
            Solver.main = Solver.numJac;
        end
    end
end