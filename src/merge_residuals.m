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
        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        % There should be a better way to do this! --> with numeric_jacobian
        if ~diff(size(J1))
            %delta = 10^-8;
            %R1p = fun(x(1:end-1),x(end)+delta);
            J1l = numeric_jacobian(@(x) fun(x,x_all(length(x)+1:end)), x, 'is', length(x));
            %J1l = (R1p-R1)*1/delta;
        else
            J1l = [];
        end
        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        varargout{1} = [R1;R2];
        varargout{2} = [J1,J1l;J2];
    else
        R1 = fun(x(1:end-1),x(end));
        R2 = res_arle(x,x_all,ds);
        varargout{1} = [R1;R2];
    end
end