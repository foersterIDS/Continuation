%% path continuation - residual_arclength_linear
%
%   Institute of Dynamics and Vibration Research
%   Leibniz University Hannover
%   17.09.2020 - Tido Kubatschek
%   22.09.2020 - Alwin Förster
%
function [residual,jacobian] = residual_arclength_linear(x,xs,ds,Opt)
    [b,a] = size(xs);    
    % approximate tangent with secant
    if a == 1
        sec = [zeros(b-1,1);1];
        sec = Opt.direction * ds * sec;
        xip1 = xs(:,end) + sec;
    else
        sec = xs(:,end) - xs(:,end-1);
        sec = ds*sec/sqrt(sum(sec.^2));
        xi = xs(:,end);
        xip1 = xi + sec;
    end
    %
    residual = sec.' * (x - xip1);
    jacobian = sec.';
end