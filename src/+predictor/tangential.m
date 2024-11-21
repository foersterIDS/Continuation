%% path continuation - predictor.tangential
%
%   Institute of Dynamics and Vibration Research
%   Leibniz University Hannover
%   12.11.2020 - Alwin Förster
%
function [funPredictor,JacPredictor] = tangential(oih,s,solverJacobian,fun,idx)
    xAll = oih.path.xAll;
    xi = xAll(:,idx);
    nd = length(xi);
    %% check jacobian:
    [nj1,nj2] = size(solverJacobian);
    if nj1>=nd-1
        if nj2==nd
            jac = solverJacobian(1:nd-1,1:nd);
        else
            jac = [solverJacobian(1:nd-1,1:nj2),aux.numericJacobian(@(x) fun(x(1:end-1),x(end)),xi,'derivativeDimensions',(nj2+1):nd,'diffquot',oih.opt.diffquot,'diffStep',oih.opt.diffStep)];
        end
    else
        jac = aux.numericJacobian(@(x) fun(x(1:end-1),x(end)),xi,'diffquot',oih.opt.diffquot,'diffStep',oih.opt.diffStep);
    end
    %% build system of equations:
    [~,m] = max(abs(diff(xAll(:,idx+(-1:0)),1,2)));
    u = 1:nd;
    u(m) = [];
    dxmds = sign(diff(xAll(m,idx+(-1:0)),1,2));
    jacM = jac(:,m);
    jacU = jac(:,u);
    % check cond.:
    if cond(jacU)<10^8
        %% tangential predictor
        dxuds = jacU\(-jacM*dxmds);
        c = dxuds/abs(dxmds);
        cex = [c(1:(m-1));dxmds;c(m:end)];
        a = sqrt(1/(sum(cex.^2)));
        dxds = cex*a;
        funPredictor = xi+dxds*s;
        JacPredictor = dxds;
    else
        %% secant predictor
        nOrder = 1;
        nFit = 0;
        if nargout>1
            [funTaylor,JacTaylor] = predictor.taylor(oih,nOrder,nFit,idx);
            JacPredictor = JacTaylor(s);
        else
            funTaylor = predictor.taylor(oih,nOrder,nFit,idx);
        end
        funPredictor = funTaylor(s);
        aux.printLine(oih,'------> Switched to basic secant predictor due to bad cond(J).\n');
    end
end