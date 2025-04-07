%% path continuation - continuation.predictor
%
%   Institute of Dynamics and Vibration Research
%   Leibniz University Hannover
%   08.05.2020 - Alwin Förster
%
function [vp,lp,funPredictor,sp,ds] = predictor(oih,ds,solverJacobian,func,resCorr,idx)
    arguments
        oih 
        ds 
        solverJacobian 
        func 
        resCorr 
        idx = []
    end
    if isempty(idx)
        idx = oih.path.nAll;
    end
    %% get funPredictor:
    if oih.opt.predictor.polynomial
        if idx==1
            funPredictor = @(s) predictor.initial(oih,s);
        else
            [nt,nf] = predictor.adaptive(oih,idx);
            [fpt,Jpt] = predictor.taylor(oih,nt,nf,idx);
            funPredictor = @(s) aux.fncHndToVal(s,fpt,Jpt);
        end
    elseif oih.opt.predictor.tangential
        if idx==1
            funPredictor = @(s) predictor.initial(oih,s);
        else
            funPredictor = @(s) predictor.tangential(oih,s,solverJacobian,func,idx);
        end
    else
        error('predictor not set or of unknown type');
    end
    %% predictorSolver:
    if oih.opt.predictorSolver
        xi = oih.path.xAll(:,idx);
        if oih.opt.enforceDsMax
            %% enforceDsMax:
            p = 0.95;
            funPred = @(s,ds) aux.mergeArlePred(funPredictor,resCorr,s,xi,ds,solverJacobian);
            funDsMax = @(s,ds) heaviside(-min(oih.opt.dsMax*p-(funPredictor(s)-xi)))*min(oih.opt.dsMax*p-(funPredictor(s)-xi));
            funSolve = @(spds) [funPred(spds(1),spds(2));funDsMax(spds(1),spds(2))];
            [spds,~,exitflag] = oih.solver.numJac(funSolve,[ds;ds]);
            if exitflag<=0
                sp = ds;
            else
                sp = spds(1);
                ds = spds(2);
            end
        else
            %% solve arc-length
            funSolve = @(s) aux.mergeArlePred(funPredictor,resCorr,s,xi,ds,solverJacobian);
            [sp,~,exitflag] = oih.solver.predictor(funSolve,ds);
            if exitflag<=0
                sp = ds;
            end
        end
    else
        sp = ds;
    end
    %% get predictor:
    xip1 = funPredictor(sp);
    %% correctPredictor:
    if oih.opt.correctPredictor && idx>1
        dxi = oih.path.xAll(:,idx)-oih.path.xAll(:,idx-1);
        dxip1 = xip1-oih.path.xAll(:,idx);
        if dot(dxip1,dxi)<0 && ~sum(ds<0)
            dxip1 = -dxip1;
            xip1 = oih.path.xAll(:,idx)+dxip1;
        end
    end
    %% make output:
    vp = xip1(1:end-1);
    lp = xip1(end);
end