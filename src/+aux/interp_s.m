%% path continuation - aux.interp_s
%
%   Institute of Dynamics and Vibration Research
%   Leibniz University Hannover
%   27.10.2020 - Tido Kubatschek
%
function [s_all_i] = interp_s(s_all, var_all, l_all, Opt)
%     var_all = vs;
%     l_all = ls;
%     s_all = 0;
%     for k = 2:length(x_all)
%         s_all = [s_all,s_all(end) + norm(x_all(:,k) - x_all(:,k-1))];
%     end

    n = length(l_all);
    
    if n < 2
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
        
    elseif Opt.interpolation_method.beziere
        aux.print_line(Opt,'--> method has not been implemented yet! Using spline instead\n');
        pp = @(s,x) spline(s, x);
    else
        error('There is no such interpolation method.')
    end
    
    % calc new s with finer inc
    
    n_points = max([n, 4]);
    
    s_i= linspace(s_all(end-n_points), s_all(end), Opt.interp_s_inc);
    x_all_i = ppval(pp(s_all, x_all), s_i);
    
    s_all_n = 0;
    
    for k = 2:length(x_all_i)
        s_all_n = [s_all_n,s_all_n(end) + norm(x_all_i(:,k) - x_all_i(:,k-1))];
    end
    
end