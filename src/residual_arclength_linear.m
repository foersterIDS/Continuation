%% path continuation - residual_arclength_linear
%
%   Institute of Dynamics and Vibration Research
%   Leibniz University Hannover
%   17.09.2020 - Tido Kubatschek
%
function [residual,jacobian] = residual_arclength_linear(x,xs,values)
    % Wichtig: Bifurkationen werden nicht richtig erkannt!!
    [b,a] = size(xs);
    
    % get information from values
    ds = values(1);
    size_l = values(2);
    
    % approximate tangent with secant
    if a == 1
        sec = [zeros(b-size_l,1);ones(size_l,1)];
        sec = ds * sec;
        xip1 = xs(:,end) + sec;
    else
        sec = xs(:,end) - xs(:,end-1);
        xi = xs(:,end);
        xip1 = xi + sec;
    end
    %
    residual = sec.' * (x - xip1);
    jacobian = sec.';
end