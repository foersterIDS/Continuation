%% path continuation - homotopy.f2
%
%   Institute of Dynamics and Vibration Research
%   Leibniz University Hannover
%   04.06.2020 - Alwin Förster
%
function [x0,l0] = f2LMInit(R,x0)
    n = length(x0);
    Mx = [ones(n-1,1),-eye(n-1)];
    Ml = [1,zeros(1,n-1)];
    Rlm = @(x) Mx*R(x).^2;
    i = 1;
    iMax = 10;
    opt = optimoptions('fsolve','Algorithm','Levenberg-Marquardt','Display','None');
    while i<=iMax
        x00 = rand*ones(n,1);
        [xs,~,exfl] = fsolve(Rlm,x00,opt);
        i = i+1;
        if i>iMax && exfl<=0
            error('No solution found');
        elseif exfl>0
            x0 = xs;
            l0 = 1+Ml*R(x0).^2;
            break;
        end
    end
end