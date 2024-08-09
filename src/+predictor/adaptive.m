%% path continuation - predictor.adaptive
%
%   Institute of Dynamics and Vibration Research
%   Leibniz University Hannover
%   18.10.2020 - Alwin Förster
%
function [nt,nf] = adaptive(oih)
    if oih.opt.predictorPolynomialAdaptive
        varOld = oih.path.varAll(:,1:end-1);
        lOld = oih.path.lAll(1:end-1);
        sOld = oih.path.sAll(1:end-1);
        varSolution = oih.path.varAll(:,end);
        lSolution = oih.path.lAll(end);
        xSolution = [varSolution;lSolution];
        errmin = inf;
        dsOld = norm(xSolution-[varOld(:,end);lOld(end)]);
        for kt=1:oih.opt.predictorPolynomialDegree
            for kf=0:oih.opt.predictorPolynomialFit
                if length(lOld)==1
                    if numel(oih.opt.direction)==1
                        xPredictorOld = [varOld;lOld+sign(oih.opt.direction)*dsOld];
                    else
                        xPredictorOld = [varOld;lOld]+oih.opt.direction*dsOld;
                    end
                else
                    oih.path.toggleStepback();
                    fpt = predictor.taylor(oih,kt,kf);
                    xPredictorOld = fpt(dsOld);
                    oih.path.toggleStepback();
                end
                err = norm(xPredictorOld-xSolution);
                if err<errmin
                    errmin = err;
                    nt = kt;
                    nf = kf;
                end
            end
        end
    else
        nt = oih.opt.predictorPolynomialDegree;
        nf = oih.opt.predictorPolynomialFit;
    end
end