function [vAll,lAll] = refinePath(fun,vAll,lAll,nRef)
    %% arguments
    arguments
        fun (1,1) function_handle
        vAll (:,:) double
        lAll (1,:) double
        nRef (1,1) double {mustBeInteger,mustBeGreaterThan(nRef,0)}
    end
    %% run
    solverOpt = optimoptions('fsolve','Display','off');
    for ii=1:nRef
        %% get idxMax
        xAll = [vAll;lAll];
        sAll = [0,cumsum(sqrt(sum(diff(xAll,1,2).^2,1)))];
        dv = sqrt(sum(diff(vAll,1,2).^2,1));
        ds = diff(sAll);
        [~,idxMax] = max(dv.*ds);
        %% predict
        sPredictor = mean(sAll(idxMax+[0,1]));
        xPredictor = spline(sAll,xAll,sPredictor);
        %% solve
        [xSolution,~,exitflag] = fsolve(@(x) [fun(x(1:(end-1)),x(end));norm(x-xAll(:,idxMax))-norm(x-xAll(:,idxMax+1))],xPredictor,solverOpt);
        dsMax = diff(sAll(idxMax+[0,1]));
        dsL = norm(xSolution-xAll(:,idxMax));
        dsU = norm(xAll(:,idxMax+1)-xSolution);
        if exitflag>0 && dsL<dsMax && dsU<dsMax
            vAll = [vAll(:,1:idxMax),xSolution(1:(end-1)),vAll(:,(idxMax+1):end)];
            lAll = [lAll(1:idxMax),xSolution(end),lAll((idxMax+1):end)];
        else
            error('unhandled error');
        end
        fprintf('refinePath: %.2f %%\n',ii/nRef*100);
    end
end