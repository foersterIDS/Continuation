%% path continuation - aux.pathFusion
%  Fuses different paths. 
% 
% 
%   Inputs:
%       vData    -- contains a vector corresponding values or a 
%                    matrix which rows are vectors of corresponding values 
%       lData    -- contains parameter points
% 
%   Outputs:
%       vFusion     -- fused vData
%       lFusion     -- fused lData
%       sFusion     -- fused sData
% 
%   Institute of Dynamics and Vibration Research
%   Leibniz University Hannover
%   28.06.2023 - Alwin Förster
%
function [vFusion,lFusion,sFusion] = pathFusion(vData,lData)
    %% arguments
    arguments (Repeating)
        vData (:,:) double
        lData (1,:) double
    end
    %% settings
    samePointTol = 10^-8;
    %% init.
    nPaths = numel(lData);
    nv = numel(vData{1}(:,1));
    nx = nv+1;
    xData = cell(1,nPaths);
    lDataStart = NaN(nPaths,1);
    for ii=1:nPaths
        if lData{ii}(end)<lData{ii}(1)
            vData{ii} = vData{ii}(:,end:-1:1);
            lData{ii} = lData{ii}(end:-1:1);
        end
        lDataStart(ii) = lData{ii}(1);
        xData{ii} = [vData{ii};lData{ii}];
    end
    [lDataStart,idxSlDS] = sort(lDataStart);
    xData = xData(idxSlDS);
    vData = vData(idxSlDS);
    lData = lData(idxSlDS);
    xFusion = xData(1);
    %% fusion
    if nPaths>1
        cFusionPath = 1; % counter for current fusion path
        openPaths = 2:nPaths;
        while ~isempty(openPaths)
            cFusionPointOfInterest = 1;
            cPathOfInterest = openPaths(1);
            xPOI = xData{cPathOfInterest};
            isMatching = false;
            doCompare = true;
            while doCompare
                %% alphaTol
                dxFusion = diff(xFusion{cFusionPath},1,2);
                alphaTol = 0;
                for ii=1:(numel(dxFusion(1,:))-1)
                    alphaTol = max(alphaTol,acos((dxFusion(:,ii+1)'*dxFusion(:,ii))/(sqrt(dxFusion(:,ii+1)'*dxFusion(:,ii+1))*sqrt(dxFusion(:,ii)'*dxFusion(:,ii)))));
                end
                alphaTol = alphaTol*1.5;
                %% current fusion path
                dx0 = xFusion{cFusionPath}(:,cFusionPointOfInterest+1)-xFusion{cFusionPath}(:,cFusionPointOfInterest);
                dsFusion = norm(dx0);
                if cFusionPointOfInterest>1
                    dx0 = (dx0+(xFusion{cFusionPath}(:,cFusionPointOfInterest)-xFusion{cFusionPath}(:,cFusionPointOfInterest-1)))/2;
                end
                %% matches?
                match = sqrt(sum((xPOI(:,[1,end])-xFusion{cFusionPath}(:,cFusionPointOfInterest)).^2,1))<dsFusion;
                if ~match(1) && match(2)
                    xPOI = xPOI(:,end:-1:1);
                    match = match(2:-1:1);
                end
                if match(1)
                    if numel(xPOI(1,:))>1
                        dxPOI = xPOI(:,2)-xPOI(:,1);
                    else
                        dxPOI = xPOI(:,1)-xFusion{cFusionPath}(:,cFusionPointOfInterest);
                    end
                    alpha = acos((dxPOI'*dx0)/(sqrt(dxPOI'*dxPOI)*sqrt(dx0'*dx0)));
                end
                if sum(match) && alpha<alphaTol
                    isMatching = true;
                end
                %% matches!
                if isMatching
                    %% fusion
                    dxPOI = xPOI(:,1)-xFusion{cFusionPath}(:,cFusionPointOfInterest);
                    sFusionTemp = [0,cumsum(sqrt(sum((xFusion{cFusionPath}(:,2:end)-xFusion{cFusionPath}(:,1:(end-1))).^2)))];
                    sPOITemp = [0,cumsum(sqrt(sum((xPOI(:,2:end)-xPOI(:,1:(end-1))).^2)))]+norm(dxPOI)+sFusionTemp(cFusionPointOfInterest);
                    xFusion{cFusionPath} = [xFusion{cFusionPath},xPOI];
                    sFusion = [sFusionTemp,sPOITemp];
                    [sFusion,idxS] = sort(sFusion);
                    xFusion{cFusionPath} = xFusion{cFusionPath}(:,idxS);
                    idxSamePoint = find(sqrt(sum(diff(xFusion{cFusionPath},1,2).^2))<samePointTol)+1;
                    if ~isempty(idxSamePoint)
                        xFusion{cFusionPath}(:,idxSamePoint) = [];
                    end
                    %% prepare for next comparison
                    if ~isempty(openPaths)
                        openPaths(1) = [];
                    end
                    cFusionPointOfInterest = 1;
                    if ~isempty(openPaths)
                        cPathOfInterest = openPaths(1);
                        xPOI = xData{cPathOfInterest};
                        isMatching = false;
                    end
                else
                    cFusionPointOfInterest = cFusionPointOfInterest+1;
                end
                %% doesn't match
                if cFusionPointOfInterest == numel(xFusion{cFusionPath}(1,:)) || isempty(openPaths)
                    doCompare = false;
                end
            end
            %% prepare for a seperate fusion path
            if ~isempty(openPaths)
                cFusionPath = cFusionPath+1;
                xFusion{cFusionPath} = xData{openPaths(1)};
                openPaths(1) = [];
            end
        end
    end
    %% write fused paths in one path
    vFusion = [];
    lFusion = [];
    sFusion = [];
    for ii=1:numel(xFusion)
        if ii==1
            sepLength = 0;
        else
            sepLength = 1;
        end
        vFusion = [vFusion,NaN(nv,sepLength),xFusion{ii}(1:nv,:)];
        lFusion = [lFusion,NaN(1,sepLength),xFusion{ii}(nx,:)];
        sFusionTemp = [0,cumsum(sqrt(sum((xFusion{ii}(:,2:end)-xFusion{ii}(:,1:(end-1))).^2)))];
        sFusion = [sFusion,NaN(1,sepLength),sFusionTemp];
    end
end