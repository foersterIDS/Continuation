%% path continuation - aux.includeReverse
%
%   Institute of Dynamics and Vibration Research
%   Leibniz University Hannover
%   05.10.2020 - Alwin Förster
%
function [PathOut,Jacobian] = includeReverse(xSolution,Path,Jacobian)
    xAll = [Path.varAll;Path.lAll];
    dist = sqrt(sum((xAll-kron(xSolution,ones(1,length(Path.lAll)))).^2));
    [~,i1] = min(dist);
    dist(i1) = inf;
    [~,i2] = min(dist);
    PathOut = Path;
    if abs(i1-i2)==1
        ii = sort([i1,i2]);
        PathOut.varAll = [Path.varAll(:,1:ii(1)),xSolution(1:end-1),Path.varAll(:,ii(2):end)];
        PathOut.lAll = [Path.lAll(1:ii(1)),xSolution(end),Path.lAll(ii(2):end)];
        Jacobian.all = cat(3,cat(3,Jacobian.all(1:ii(1)),Jacobian.last),Jacobian.all(:,:,ii(2):end));
    else
        PathOut.varAll = Path.varAll;
        PathOut.lAll = Path.lAll;
    end
end