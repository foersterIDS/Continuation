%% path continuation - continuation.predictor
%
%   Institute of Dynamics and Vibration Research
%   Leibniz University Hannover
%   08.05.2020 - Alwin Förster
%
function [vp,lp,funPredictor,sp,ds] = predictor(Path,ds,solverJacobian,func,resCorr,Solver,Opt)
    %% get funPredictor:
    if Opt.predictor.polynomial
        if numel(Path.lAll)==1
            funPredictor = @(s) predictor.initial(Path,s,Opt);
        else
            [nt,nf] = predictor.adaptive(Path,Opt);
            [fpt,Jpt] = predictor.taylor(Path,nt,nf);
            funPredictor = @(s) aux.fncHndToVal(s,fpt,Jpt);
        end
    elseif Opt.predictor.tangential
        if numel(Path.lAll)==1
            funPredictor = @(s) predictor.initial(Path,s,Opt);
        else
            funPredictor = @(s) predictor.ode(Path,s,solverJacobian,func,Opt);
        end
    else
        error('predictor not set or of unknown type');
    end
    %% predictorSolver:
    if Opt.predictorSolver
        xi = [Path.varAll(:,end);Path.lAll(end)];
        if Opt.enforceDsMax
            %% enforceDsMax:
            p = 0.95;
            funPred = @(s,ds) aux.mergeArlePred(funPredictor,resCorr,s,xi,ds,solverJacobian);
            funDsMax = @(s,ds) heaviside(-min(Opt.dsMax*p-(funPredictor(s)-xi)))*min(Opt.dsMax*p-(funPredictor(s)-xi));
            funSolve = @(spds) [funPred(spds(1),spds(2));funDsMax(spds(1),spds(2))];
            [spds,~,exitflag] = Solver.numJac(funSolve,[ds;ds]);
            if exitflag<=0
                sp = ds;
            else
                sp = spds(1);
                ds = spds(2);
            end
        else
            %% solve arc-length
            funSolve = @(s) aux.mergeArlePred(funPredictor,resCorr,s,xi,ds,solverJacobian);
            [sp,~,exitflag] = Solver.predictor(funSolve,ds);
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
    if Opt.correctPredictor && numel(Path.lAll)>1
        dxi = [Path.varAll(:,end);Path.lAll(end)]-[Path.varAll(:,end-1);Path.lAll(end-1)];
        dxip1 = xip1-[Path.varAll(:,end);Path.lAll(end)];
        if dot(dxip1,dxi)<0 && ~sum(ds<0)
            dxip1 = -dxip1;
            xip1 = [Path.varAll(:,end);Path.lAll(end)]+dxip1;
        end
    end
    %% make output:
    vp = xip1(1:end-1);
    lp = xip1(end);
end