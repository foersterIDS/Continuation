%% path continuation - residual_corrector_parabola
%
%   Institute of Dynamics and Vibration Research
%   Leibniz University Hannover
%   29.12.2021 - Alwin Förster
%
function [residual,jacobian] = residual_corrector_parabola(x,x_all,ds,Opt)
    [b,a] = size(x_all);    
    %% approximate tangent with secant
    %
    if a == 1
        if numel(Opt.direction)==1
            sec = [zeros(b-1,1);1];
            sec = Opt.direction * ds * sec;
        else
            sec = Opt.direction * ds;
        end
        xip1 = x_all(:,end) + sec;
    else
        sec = x_all(:,end) - x_all(:,end-1);
        sec = ds*sec/sqrt(sum(sec.^2));
        xi = x_all(:,end);
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
    r_para = -sum((chi-xip1).^2./sf);
    x_para = chi+sec*r_para;
    %
    %% calc. residual
    %
    residual = norm(x_para-x);
    jacobian = norm(sec)*(-(sec.'/(sec.'*sec))-sum((2*diag(x+sec*(-(sec.'/(sec.'*sec))*(x-xip1))-xip1)*(eye(b)+sec*(-(sec.'/(sec.'*sec))*eye(b))))./sf));
    %
end