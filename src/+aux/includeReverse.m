%% path continuation - aux.includeReverse
%
%   Institute of Dynamics and Vibration Research
%   Leibniz University Hannover
%   05.10.2020 - Alwin Förster
%
function includeReverse(xSolution,oih)
    xAll = oih.path.xAll;
    dist = sqrt(sum((xAll-kron(xSolution,ones(1,oih.path.nAll))).^2));
    [~,i1] = min(dist);
    dist(i1) = inf;
    [~,i2] = min(dist);
    if abs(i1-i2)==1
        ii = sort([i1,i2]);
        oih.path.addPoint(xSolution(1:end-1),xSolution(end),oih.solver.jacobian);
    end
end