%% path continuation - corrector.residualParaboloid
%
%   Institute of Dynamics and Vibration Research
%   Leibniz University Hannover
%   29.12.2021 - Alwin Förster
%
function [residual,jacobian] = residualParaboloid(x,xAll,ds,Opt)
    [b,a] = size(xAll);    
    %% approximate tangent with secant
    %
    if a == 1
        if numel(Opt.direction)==1
            sec = [zeros(b-1,1);1];
            sec = Opt.direction * ds * sec;
        else
            sec = Opt.direction * ds;
        end
        xip1 = xAll(:,end) + sec;
    else
        sec = xAll(:,end) - xAll(:,end-1);
        sec = ds*sec/sqrt(sum(sec.^2));
        xi = xAll(:,end);
        xip1 = xi + sec;
    end
    %
    %% projection to orthogonal plane
    %
    r = -(sec.'*(x-xip1))/(sec.'*sec);
    chi = x+sec*r; % projection of x to orthogonal plane
    %
    %% calc. squared value
    %
    sf = (ds/10)^2/ds; % scaling factor
    rPara = -sum((chi-xip1).^2./sf);
    xPara = chi+sec*rPara;
    %
    %% calc. residual
    %
    residual = norm(xPara-x);
    jacobian = norm(sec)*(-(sec.'/(sec.'*sec))-sum((2*diag(x+sec*(-(sec.'/(sec.'*sec))*(x-xip1))-xip1)*(eye(b)+sec*(-(sec.'/(sec.'*sec))*eye(b))))./sf));
    %
end