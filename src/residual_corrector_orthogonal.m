%% path continuation - residual_corrector_orthogonal
%
%   Institute of Dynamics and Vibration Research
%   Leibniz University Hannover
%   17.09.2020 - Tido Kubatschek
%   22.09.2020 - Alwin Förster
%
function [residual,jacobian] = residual_corrector_orthogonal(x,x_all,ds,Opt)
    [b,a] = size(x_all);    
    % approximate tangent with secant
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
    residual = sec.' * (x - xip1);
    jacobian = sec.';
end