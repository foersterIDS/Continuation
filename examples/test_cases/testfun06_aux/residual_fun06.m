% function [f,J] = residual_fun06(v,l,Nres)
%     f = [v(1)^2-1;v(2:end)+l*v(1:end-1).^5];
%     J = [2*v(1),zeros(1,Nres-1);
%          [zeros(Nres-1,1),eye(Nres-1)]+5*l*[diag(v(1:end-1).^4),zeros(Nres-1,1)]];
% end

function [f,J] = residual_fun06(v,l)
    f = [v(1)-l.^2; v(2) - 10^-4 + 2*10^-5*sin(l)];
    J = [1, 0, 2*l; 
        0, 1, 2*10^-5*cos(l)];
end
