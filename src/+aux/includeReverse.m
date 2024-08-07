%% path continuation - aux.includeReverse
%
%   Institute of Dynamics and Vibration Research
%   Leibniz University Hannover
%   05.10.2020 - Alwin Förster
%
function [Path,Jacobian] = includeReverse(xSolution,Path,Jacobian)
    xAll = [Path.varAll;Path.lAll];
    dist = sqrt(sum((xAll-kron(xSolution,ones(1,length(Path.lAll)))).^2));
    [~,i1] = min(dist);
    dist(i1) = inf;
    [~,i2] = min(dist);
    if abs(i1-i2)==1
        ii = sort([i1,i2]);
        Path.addPoint(xSolution(1:end-1),xSolution(end),Jacobian.last);
        Jacobian.all = cat(3,cat(3,Jacobian.all(1:ii(1)),Jacobian.last),Jacobian.all(:,:,ii(2):end));
    end
end