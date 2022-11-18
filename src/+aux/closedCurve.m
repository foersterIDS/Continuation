%% path continuation - aux.closedCurve
%
%   Institute of Dynamics and Vibration Research
%   Leibniz University Hannover
%   03.11.2020 - Tido Kubatschek
%   09.11.2020 - Alwin Förster
%
function [isClosed, Opt, Counter] = closedCurve(Opt,Path,ds,Counter)
    isClosed = 0;
    n = 5;
    np = 4; % polynomial order (np = 2*k with k in N)
    epsPol = 2.5e-2;
    nu = 1+np/2;
    ns = numel(Path.lAll);
    if ns > n+nu
        xAll = [Path.varAll; Path.lAll];
        epsDistX = min([2*norm(xAll(:,end)-xAll(:,end-1)),2*ds,Opt.dsMax]);
        epsDistL = min([2*norm(xAll(end,end)-xAll(end,end-1)),2*ds,Opt.dsMax]);
        epsAng = 0.5*(2*pi / 360);
        epsDir = 1e-3;
        %% exclude points before current points which are too close
        %
        ignoredDist = 2*norm(xAll(:,end)-xAll(:,end-1));
        distanceX = sqrt(sum((xAll(:,1:(ns - 1)) - xAll(:,end)).^2,1));
        distanceL = sqrt((xAll(end,1:(ns - 1)) - xAll(end,end)).^2);
        lastInd = find(distanceX > ignoredDist);
        if isempty(lastInd)
            flag = 0;
        else
            ignored = ns - lastInd(end);
            %
            flag = 0;
            if ns>n
                distX = distanceX(nu:(ns - ignored));
                distL = distanceL(nu:(ns - ignored));
                kFlags = find((distX <= epsDistX).*(distL <= epsDistL));
                if numel(kFlags)>0
                    flag = 1;
                    kFlags = kFlags + nu - 1;
                end
            end
        end       
        %% check wether the found point is crossed under the same angle
        countFlag = 0;
        if flag == 1
            for i=1:numel(kFlags)
                %% polynomials
                % centered polynomial to matching point
                s0F = Path.sAll(kFlags(i));
                pF = poly.fitn(Path.sAll(kFlags(i)+(-np/2:np/2))-s0F,xAll(:,kFlags(i)+(-np/2:np/2)),np);
                % arclength correction
                dsc = distX(kFlags(i)-nu+1);
                % centered polynomials to last point
                s0C = Path.sAll(end);
                pCP = poly.fitn(Path.sAll(end+(-np:0))-s0C+dsc,xAll(:,end+(-np:0)),np);
                pCM = poly.fitn(Path.sAll(end+(-np:0))-s0C-dsc,xAll(:,end+(-np:0)),np);
                % check polyinomials
                if min([norm(pF-pCP),norm(pF-pCM)]/norm(pF)) <= epsPol
                    countFlag = 1;
                    break;
                end
                %% angle
                vecF = xAll(:,kFlags(i)) - xAll(:,kFlags(i)-1);
                vecC = xAll(:,end) - xAll(:,end - 1);
                angle = aux.vectorAngle(vecF, vecC);
                %% direction
                rF = vecF / norm(vecF);
                rC = vecC / norm(vecC);
                dir = norm(rF - rC);
                %% check
                if angle <= epsAng && dir <= epsDir % check angle and direction
                    countFlag = 1;
                end
            end
        end
        if countFlag
            Counter.closedCurve = Counter.closedCurve + 1;
        elseif Counter.closedCurve > 0
            Counter.closedCurve = Counter.closedCurve - 1;
        end
        if Counter.closedCurve == Opt.maxClosedCounter
            isClosed = 1;
        end
    end
end

