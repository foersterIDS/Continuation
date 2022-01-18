%% path continuation - predictor.ode
%
%   Institute of Dynamics and Vibration Research
%   Leibniz University Hannover
%   12.11.2020 - Alwin Förster
%
function [fun_predictor,Jac_predictor] = ode(Path,s,solver_jacobian,fun)
    x_all = [Path.var_all;Path.l_all];
    xi = x_all(:,end);
    nd = length(xi);
    %% check jacobian:
    [nj1,nj2] = size(solver_jacobian);
    if nj1>=nd-1
        if nj2==nd
            jac = solver_jacobian(1:nd-1,1:nd);
        else
            jac = [solver_jacobian(1:nd-1,1:nj2),numeric_jacobian(@(x) fun(x(1:end-1),x(end)),xi,'derivative_dimensions',(nj2+1):nd,'diffquot',Opt.diffquot)];
        end
    else
        jac = numeric_jacobian(@(x) fun(x(1:end-1),x(end)),xi,'diffquot',Opt.diffquot);
    end
    %% build system of equations:
    [~,m] = max(abs(diff(x_all(:,end+(-1:0)),1,2)));
    u = 1:nd;
    u(m) = [];
    dxmds = sign(diff(x_all(m,end+(-1:0)),1,2));
    jac_m = jac(:,m);
    jac_u = jac(:,u);
    dxuds = jac_u\(-jac_m*dxmds);
    c = dxuds/abs(dxmds);
    cex = [c(1:(m-1));dxmds;c(m:end)];
    a = sqrt(1/(sum(cex.^2)));
    dxds = cex*a;
    fun_predictor = xi+dxds*s;
    Jac_predictor = dxds;
end