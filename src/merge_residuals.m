%% path continuation - merge_residuals
%
%   Institute of Dynamics and Vibration Research
%   Leibniz University Hannover
%   19.05.2020 - Alwin Förster
%   16.09.2020 - Tido Kubatschek 
%
function [varargout] = merge_residuals(fun,res_arle,x,x_all,ds,Opt)
    if Opt.jacobian
        [R1,J1] = fun(x(1:end-1),x(end));
        [R2,J2] = res_arle(x,x_all,ds);
        if ~diff(size(J1))
            J1l = numeric_jacobian(@(x) fun(x(1:end-1),x(end)), x, 'derivative_dimensions', length(x), 'jacob', Opt.jacob);
        else
            J1l = [];
        end
        varargout{1} = [R1;R2];
        varargout{2} = [J1,J1l;J2];
    else
        R1 = fun(x(1:end-1),x(end));
        R2 = res_arle(x,x_all,ds);
        varargout{1} = [R1;R2];
    end
end