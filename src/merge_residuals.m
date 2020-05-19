%% path continuation - merge_residuals
%
%   Institute of Dynamics and Vibration Research
%   Leibniz University Hannover
%   19.05.2020 - Alwin Förster
%
function [varargout] = merge_residuals(fun,res_arle,x,x_all,ds)
    is_jacobian_1 = false;
    is_jacobian_2 = false;
    if abs(nargout(fun))==2
        [R1,J1] = fun(x(1:end-1),x(end));
        is_jacobian_1 = true;
    else
        R1 = fun(x(1:end-1),x(end));
    end
    if abs(nargout(res_arle))==2
        [R2,J2] = res_arle(x,x_all,ds);
        is_jacobian_2 = true;
    else
        R2 = res_arle(x,x_all,ds);
    end
    if is_jacobian_1 && is_jacobian_2
        varargout{1} = [R1;R2];
        varargout{2} = [J1;J2];
    else
        varargout{1} = [R1;R2];
    end
end