%% path continuation - predictor.taylor
%
%   Institute of Dynamics and Vibration Research
%   Leibniz University Hannover
%   05.10.2020 - Alwin Förster
%
function [funPredictor,JacPredictor] = taylor(oih,nOrder,nFit,idx)
    nOrder = min([idx-1,nOrder]);
    ns = min([idx,nOrder+1+nFit]);
    xAll = oih.path.xAll;
    nd = length(xAll(:,idx));
    %% calc scaling:
    dscMin = 10^-15;
    xBasis = xAll(:,idx-ns+1);
    sBasis = oih.path.sAll(idx-ns+1);
    if idx>2
        dscX = max([mean(diff(abs(xAll(:,idx+((-ns+1):0))-xBasis),1,2),2),ones(nd,1)*dscMin]')';
        dscS = max([mean(diff(abs(oih.path.sAll(idx+((-ns+1):0))-sBasis),1,2),2),dscMin]);
    else
        dscX = max([abs(xAll(:,2)-xBasis),ones(nd,1)*dscMin]')';
        dscS = max([abs(oih.path.sAll(2)-sBasis),dscMin]);
    end
    %% calc taylor-predictor:
    pSc = poly.fitn((oih.path.sAll(idx+((-ns+1):0))-sBasis)./dscS,(xAll(:,idx+((-ns+1):0))-xBasis)./kron(dscX,ones(1,ns)),nOrder);
    funPredictorSc = @(s) poly.valn(pSc,(oih.path.sAll(idx)+s-sBasis)/dscS,nd);
    if nargout>1
        dpSc = poly.dern(pSc,nd);
        JacPredictorSc = @(s) poly.valn(dpSc,(oih.path.sAll(idx)+s-sBasis)/dscS,nd);
    end
    %% descale predictor:
    funPredictor = @(s) funPredictorSc(s).*dscX+xBasis;
    if nargout>1
        JacPredictor = @(s) JacPredictorSc(s).*dscX;
    end
end