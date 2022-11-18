%% path continuation - predictor.ode
%
%   Institute of Dynamics and Vibration Research
%   Leibniz University Hannover
%   12.11.2020 - Alwin Förster
%
function [funPredictor,JacPredictor] = ode(Path,s,solverJacobian,fun,Opt)
    xAll = [Path.varAll;Path.lAll];
    xi = xAll(:,end);
    nd = length(xi);
    %% check jacobian:
    [nj1,nj2] = size(solverJacobian);
    if nj1>=nd-1
        if nj2==nd
            jac = solverJacobian(1:nd-1,1:nd);
        else
            jac = [solverJacobian(1:nd-1,1:nj2),numericJacobian(@(x) fun(x(1:end-1),x(end)),xi,'derivativeDimensions',(nj2+1):nd,'diffquot',Opt.diffquot)];
        end
    else
        jac = aux.numericJacobian(@(x) fun(x(1:end-1),x(end)),xi,'diffquot',Opt.diffquot);
    end
    %% build system of equations:
    [~,m] = max(abs(diff(xAll(:,end+(-1:0)),1,2)));
    u = 1:nd;
    u(m) = [];
    dxmds = sign(diff(xAll(m,end+(-1:0)),1,2));
    jacM = jac(:,m);
    jacU = jac(:,u);
    dxuds = jacU\(-jacM*dxmds);
    c = dxuds/abs(dxmds);
    cex = [c(1:(m-1));dxmds;c(m:end)];
    a = sqrt(1/(sum(cex.^2)));
    dxds = cex*a;
    funPredictor = xi+dxds*s;
    JacPredictor = dxds;
end