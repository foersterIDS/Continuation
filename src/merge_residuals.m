%% path continuation - merge_residuals
%
%   Institute of Dynamics and Vibration Research
%   Leibniz University Hannover
%   19.05.2020 - Alwin Förster
%   16.09.2020 - Tido Kubatschek 
%
function [varargout] = merge_residuals(fun,res_corr,x,x_all,ds,last_jacobian,Opt)
    if Opt.jacobian
        [R1,J1] = fun(x(1:end-1),x(end));
        if Opt.corrector.orthogonal
            [R2,J2] = res_corr(x,x_all,ds,last_jacobian);
        else
            [R2,J2] = res_corr(x,x_all,ds);
        end
        
        [n11,n12] = size(J1);
        if n12<=n11
            J1l = numeric_jacobian(@(x) fun(x(1:end-1),x(end)), x, 'derivative_dimensions', (n12+1):numel(x), 'diffquot', Opt.diffquot, 'central_value', R1);
        else
            J1l = [];
        end
        varargout{1} = [R1;R2];
        varargout{2} = [J1,J1l;J2];
    else
        R1 = fun(x(1:end-1),x(end));
        if Opt.corrector.orthogonal
            R2 = res_corr(x,x_all,ds,last_jacobian);
        else
            R2 = res_corr(x,x_all,ds);
        end
        varargout{1} = [R1;R2];
    end
end