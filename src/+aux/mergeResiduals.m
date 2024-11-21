%% path continuation - aux.mergeResiduals
%
%   Institute of Dynamics and Vibration Research
%   Leibniz University Hannover
%   19.05.2020 - Alwin Förster
%   16.09.2020 - Tido Kubatschek 
%
function [varargout] = mergeResiduals(func,resCorr,x,xAll,ds,lastJacobian,oih)
    %% eval and merge functions
    %
    if oih.opt.jacobian
        [R1,J1] = func(x(1:end-1),x(end));
        [R2,J2] = resCorr(x,xAll,ds,lastJacobian);
        [n11,n12] = size(J1);
        [n21,n22] = size(J2);
        if n12<=n11
            J1l = aux.numericJacobian(@(x) func(x(1:end-1),x(end)), x, 'derivativeDimensions', (n12+1):numel(x), 'diffquot', oih.opt.diffquot, 'centralValue', R1,'diffStep',oih.opt.diffStep);
        elseif n12>n22
            J1 = J1(:,1:n22);
            J1l = [];
        else
            J1l = [];
        end
        R = [R1;R2];
        varargout{1} = R;
        J = [J1,J1l;J2];
        varargout{2} = J;
    else
        R1 = func(x(1:end-1),x(end));
        R2 = resCorr(x,xAll,ds,lastJacobian);
        R = [R1;R2];
        varargout{1} = R;
    end
    %
    %% preconditioning
    %
    if aux.ison(oih.opt.preconditioning)
        nR = numel(varargout{1});
        A = eye(nR);
        if oih.do.precon
            if oih.opt.preconditioning.jacobian
                if oih.opt.jacobian
                    if cond(J(1:nR,1:nR))<10^8
                        A = J(1:nR,1:nR);
                        invA = A\eye(nR);
                        varargout{1} = invA*R;
                        varargout{2} = invA*J;
                    else
                        1;
                    end
                else
                    [~,method] = aux.ison(oih.opt.preconditioning);
                    error("'jacobian' has to be 'on' to use %s preconditioning.",method);
                end
            elseif oih.opt.preconditioning.jacobiScaled
                if oih.opt.jacobian
                    e = 3;
                    a = sqrt(sum(abs(J),2).^2);
                    a = min(a,10^+e);
                    a = max(a,10^-e);
                    A = diag(a);
                    invA = diag(1./a);
                    varargout{1} = invA*R;
                    varargout{2} = invA*J;
                else
                    [~,method] = aux.ison(oih.opt.preconditioning);
                    error("'jacobian' has to be 'on' to use %s preconditioning.",method);
                end
            elseif oih.opt.preconditioning.diagJacobian
                if oih.opt.jacobian
                    e = 8;
                    a = diag(J);
                    a = min(a,10^+e);
                    a = max(a,10^-e);
                    A = diag(a);
                    invA = diag(1./a);
                    varargout{1} = invA*R;
                    varargout{2} = invA*J;
                else
                    [~,method] = aux.ison(oih.opt.preconditioning);
                    error("'jacobian' has to be 'on' to use %s preconditioning.",method);
                end
            elseif oih.opt.preconditioning.incompleteLU
                if oih.opt.jacobian
                    [L,U,~] = lu(J(:,1:nR));
                    A = L*U;
                    while cond(A)>10^8
                        A = A+eye(nR);
                    end
                    invA = A\eye(nR);
                    varargout{1} = invA*R;
                    varargout{2} = invA*J;
                else
                    [~,method] = aux.ison(oih.opt.preconditioning);
                    error("'jacobian' has to be 'on' to use %s preconditioning.",method);
                end
            else
                [~,preconMethod] = aux.ison(oih.opt.preconditioning);
                error('Preconditioning method %s not implemented!',preconMethod);
            end
        end
        oih.temp.invPrecon = A;
    end
    %
end