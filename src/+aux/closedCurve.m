%% path continuation - aux.closedCurve
%
%   Institute of Dynamics and Vibration Research
%   Leibniz University Hannover
%   03.11.2020 - Tido Kubatschek
%   09.11.2020 - Alwin Förster
%
function [isClosed] = closedCurve(oih,ds)
    isClosed = 0;
    n = 5;
    np = 4; % polynomial order (np = 2*k with k in N)
    epsPol = 2.5e-2;
    dsTol = 1.1;
    nu = 1+np/2;
    ns = oih.path.nAll;
    if ns > n+nu
        xAll = oih.path.xAll;
        epsDistX = min([dsTol*norm(xAll(1:(end-1),end)-xAll(1:(end-1),end-1)),dsTol*ds,oih.opt.dsMax]);
        epsDistL = min([dsTol*norm(xAll(end,end)-xAll(end,end-1)),dsTol*ds,oih.opt.dsMax]);
        epsAng = 0.5*(2*pi / 360);
        epsDir = 1e-3;
        %% exclude points before current points which are too close
        %
        ignoredDist = dsTol*norm(xAll(:,end)-xAll(:,end-1));
        distanceX = sqrt(sum((xAll(1:(end-1),1:(ns - 1)) - xAll(1:(end-1),end)).^2,1));
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
                s0F = oih.path.sAll(kFlags(i));
                pF = poly.fitn(oih.path.sAll(kFlags(i)+(-np/2:np/2))-s0F,xAll(:,kFlags(i)+(-np/2:np/2)),np);
                % arclength correction
                dsc = distX(kFlags(i)-nu+1);
                % centered polynomials to last point
                s0C = oih.path.sAll(end);
                pCP = poly.fitn(oih.path.sAll(end+(-np:0))-s0C+dsc,xAll(:,end+(-np:0)),np);
                pCM = poly.fitn(oih.path.sAll(end+(-np:0))-s0C-dsc,xAll(:,end+(-np:0)),np);
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
            oih.counter.closedCurve = oih.counter.closedCurve + 1;
        elseif oih.counter.closedCurve > 0
            oih.counter.closedCurve = oih.counter.closedCurve - 1;
        end
        if oih.counter.closedCurve == oih.opt.maxClosedCounter
            isClosed = 1;
        end
    end
end

