%% path continuation - predictor.taylor
%
%   Institute of Dynamics and Vibration Research
%   Leibniz University Hannover
%   05.10.2020 - Alwin Förster
%
function [funPredictor,JacPredictor] = taylor(oih,no,nf)
    no = min([oih.path.nAll-1,no]);
    ns = min([oih.path.nAll,no+1+nf]);
    xAll = oih.path.xAll;
    nd = length(xAll(:,end));
    %% calc scaling:
    dscMin = 10^-15;
    xBasis = xAll(:,end-ns+1);
    sBasis = oih.path.sAll(end-ns+1);
    if oih.path.nAll>2
        dscX = max([mean(diff(abs(xAll(:,end+((-ns+1):0))-xBasis),1,2),2),ones(nd,1)*dscMin]')';
        dscS = max([mean(diff(abs(oih.path.sAll(end+((-ns+1):0))-sBasis),1,2),2),dscMin]);
    else
        dscX = max([abs(xAll(:,2)-xBasis),ones(nd,1)*dscMin]')';
        dscS = max([abs(oih.path.sAll(2)-sBasis),dscMin]);
    end
    %% calc taylor-predictor:
    pSc = poly.fitn((oih.path.sAll(end+((-ns+1):0))-sBasis)./dscS,(xAll(:,end+((-ns+1):0))-xBasis)./kron(dscX,ones(1,ns)),no);
    funPredictorSc = @(s) poly.valn(pSc,(oih.path.sAll(end)+s-sBasis)/dscS,nd);
    if nargout>1
        dpSc = poly.dern(pSc,nd);
        JacPredictorSc = @(s) poly.valn(dpSc,(oih.path.sAll(end)+s-sBasis)/dscS,nd);
    end
    %% descale predictor:
    funPredictor = @(s) funPredictorSc(s).*dscX+xBasis;
    if nargout>1
        JacPredictor = @(s) JacPredictorSc(s).*dscX;
    end
end