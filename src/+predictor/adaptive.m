%% path continuation - predictor.adaptive
%
%   Institute of Dynamics and Vibration Research
%   Leibniz University Hannover
%   18.10.2020 - Alwin Förster
%
function [nt,nf] = adaptive(Path,Opt)
    if Opt.predictorPolynomialAdaptive
        varOld = Path.varAll(:,1:end-1);
        lOld = Path.lAll(1:end-1);
        sOld = Path.sAll(1:end-1);
        varSolution = Path.varAll(:,end);
        lSolution = Path.lAll(end);
        xSolution = [varSolution;lSolution];
        errmin = inf;
        dsOld = norm(xSolution-[varOld(:,end);lOld(end)]);
        for kt=1:Opt.predictorPolynomialDegree
            for kf=0:Opt.predictorPolynomialFit
                if length(lOld)==1
                    if numel(Opt.direction)==1
                        xPredictorOld = [varOld;lOld+sign(Opt.direction)*dsOld];
                    else
                        xPredictorOld = [varOld;lOld]+Opt.direction*dsOld;
                    end
                else
                    fpt = predictor.taylor(struct('varAll',varOld,'lAll',lOld,'sAll',sOld),kt,kf);
                    xPredictorOld = fpt(dsOld);
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
        nt = Opt.predictorPolynomialDegree;
        nf = Opt.predictorPolynomialFit;
    end
end