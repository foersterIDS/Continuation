%% path continuation - predictor_ode
%
%   Institute of Dynamics and Vibration Research
%   Leibniz University Hannover
%   12.11.2020 - Alwin Förster
%
function [xp] = predictor_ode(var_all,l_all,ds,solver_jacobian,fun)
    x_all = [var_all;l_all];
    xi = x_all(:,end);
    nd = length(xi);
    %% check jacobian:
    [nj1,nj2] = size(solver_jacobian);
    if nj1>=nd-1
        if nj2==nd
            jac = solver_jacobian(1:nd-1,1:nd);
        else
            jac = [solver_jacobian(1:nd-1,1:nj2),numeric_jacobian(@(x) fun(x(1:end-1),x(end)),xi,'derivative_dimensions',(nj2+1):nd)];
        end
    else
        error('wrong jacobian');
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
    xp = xi+dxds*ds;
    
    
    
%     E = eye(nd);
%     for k=nd:-1:1
%         e_k = E(:,k);
%         M = [jac;e_k.'];
%         if rank(M)==nd
%             break;
%         end
%     end
%     e_nd = E(:,nd);
%     %% calc. gradient:
%     dx = M\e_nd;
%     z = dx/norm(dx);
%     %% calc. predictor:
%     xp = xi+z*ds;
end