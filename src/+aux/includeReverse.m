%% path continuation - aux.includeReverse
%
%   Institute of Dynamics and Vibration Research
%   Leibniz University Hannover
%   05.10.2020 - Alwin Förster
%
function includeReverse(xSolution,Path,Solver)
    xAll = [Path.varAll;Path.lAll];
    dist = sqrt(sum((xAll-kron(xSolution,ones(1,Path.nAll))).^2));
    [~,i1] = min(dist);
    dist(i1) = inf;
    [~,i2] = min(dist);
    if abs(i1-i2)==1
        ii = sort([i1,i2]);
        Path.addPoint(xSolution(1:end-1),xSolution(end),Solver.jacobian);
    end
end