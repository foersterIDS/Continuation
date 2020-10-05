%% path continuation - include_reverse
%
%   Institute of Dynamics and Vibration Research
%   Leibniz University Hannover
%   05.10.2020 - Alwin Förster
%
function [var_all_out,l_all_out] = include_reverse(x_solution,var_all,l_all)
    x_all = [var_all;l_all];
    dist = sqrt(sum((x_all-kron(x_solution,ones(1,length(l_all)))).^2));
    [~,i1] = min(dist);
    dist(i1) = inf;
    [~,i2] = min(dist);
    if abs(i1-i2)==1
        ii = sort([i1,i2]);
        var_all_out = [var_all(:,1:ii(1)),x_solution(1:end-1),var_all(:,ii(2):end)];
        l_all_out = [l_all(1:ii(1)),x_solution(end),l_all(ii(2):end)];
    else
        var_all_out = var_all;
        l_all_out = l_all;
    end
end