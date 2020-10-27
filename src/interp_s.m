%% path continuation - interp_s
%
%   Institute of Dynamics and Vibration Research
%   Leibniz University Hannover
%   27.10.2020 - Tido Kubatschek
%
function [var_all_i, l_all_i, s_all_i] = interp_s(s_all, var_all, l_all, Opt)
    % determine interpolation method from Opt
    method = Opt.interpolation_method;
    
    x_all = [var_all; l_all];
    
    % calc new s with finer inc
    s_inc = Opt.interp_s_inc;
    s_all_i = s_all(1):s_inc:s_all(end);
    
    if Opt.interpolation_method.spline
        % use MATLAB intern function spline
        pp = spline(s_all, x_all);

    elseif Opt.interpolation_method.pchip
        % use MATLAB intern function pchip
        pp = pchip(s_all, x_all);
        
    elseif Opt.interpolation_method.makima
        % use MATLAB intern function makima
        pp = makima(s_all, x_all);
        
    elseif Opt.interpolation_method.beziere
        warning('method has not been implemented yet! Using spline instead');
        pp = spline(s_all, x_all);
    else
        error('There is no such interpolation method.')
    end
    
    x_all_i = ppval(pp, s_all_i);
    
    var_all_i = x_all_i(1:end-1,:);
    l_all_i = x_all_i(end,:);
end

