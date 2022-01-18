%% path continuation - aux.include_reverse
%
%   Institute of Dynamics and Vibration Research
%   Leibniz University Hannover
%   05.10.2020 - Alwin Förster
%
function [Path_out] = include_reverse(x_solution,Path)
    x_all = [Path.var_all;Path.l_all];
    dist = sqrt(sum((x_all-kron(x_solution,ones(1,length(Path.l_all)))).^2));
    [~,i1] = min(dist);
    dist(i1) = inf;
    [~,i2] = min(dist);
    Path_out = Path;
    if abs(i1-i2)==1
        ii = sort([i1,i2]);
        Path_out.var_all = [Path.var_all(:,1:ii(1)),x_solution(1:end-1),Path.var_all(:,ii(2):end)];
        Path_out.l_all = [Path.l_all(1:ii(1)),x_solution(end),Path.l_all(ii(2):end)];
    else
        Path_out.var_all = Path.var_all;
        Path_out.l_all = Path.l_all;
    end
end