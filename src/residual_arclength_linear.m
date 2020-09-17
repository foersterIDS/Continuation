%% path continuation - residual_arclength_linear
%
%   Institute of Dynamics and Vibration Research
%   Leibniz University Hannover
%   17.09.2020 - Tido Kubatschek
%
function [residual,jacobian] = residual_arclength_linear(x,xs,ds,Opt)
    
    % testing for jacobian
    
    if Opt.jacobian
        % tangent
    else
        % secant
    end
    
%     residual = ds^2-(x-xs(:,end))'*(x-xs(:,end));
%     jacobian = -2*(x-xs(:,end))';

end