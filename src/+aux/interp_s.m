%% path continuation - aux.interp_s
%
%   Institute of Dynamics and Vibration Research
%   Leibniz University Hannover
%   27.10.2020 - Tido Kubatschek
%
function [s_all_i] = interp_s(s_all, var_all, l_all, Opt)
    %%
    Opt.interpolation_method.spline = true;
    Opt.number_interp_s = 100;
    number_of_points = length(l_all);
    
    if number_of_points < 2
        error('for interpolation at least 2 data points are needed.')
    end
    
    x_all = [var_all; l_all];
    
    if Opt.interpolation_method.spline
        % use MATLAB intern function spline
        pp = @(s,x) spline(s, x);

    elseif Opt.interpolation_method.pchip
        % use MATLAB intern function pchip
        pp = @(s,x) pchip(s, x);
        
    elseif Opt.interpolation_method.makima
        % use MATLAB intern function makima
        pp = @(s,x) makima(s, x);
        
    else
        error('There is no such interpolation method.')
    end
    
    % calc new s with finer inc
    
    number_of_points = min([number_of_points, 4]);
    s_old = 0;
    s_start = s_all(end-number_of_points);
    %%
    s_all_tmp = s_all;
    while true
        s_i = linspace(s_all_tmp(end-number_of_points), s_all_tmp(end), Opt.number_interp_s);
        x_all_i = ppval(pp(s_all_tmp, x_all), s_i);

        s_all_n = s_start;

        for k = 2:length(x_all_i)
            s_all_n = [s_all_n,s_all_n(end) + norm(x_all_i(:,k) - x_all_i(:,k-1))];
        end
        
        if abs(s_all_n(end) - s_old) <= 0.1
            break;
        end
        
        s_old = s_all_n(end);
        x_all = x_all(:,1:end);
        s_all_tmp = s_all_tmp(1:end);
        
        x_all = [x_all, x_all_i(:,2:end)];
        s_all_tmp = [s_all_tmp, s_all_n(2:end)];
    end
    s_all_tmp = s_all_tmp(1:end);
    s_all_tmp = [s_all_tmp, s_all_n(2:end)];
    s_all_i = s_all_tmp;
    
    %%
end