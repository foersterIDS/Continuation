%% path continuation - aux.exitLoop
%
%   Institute of Dynamics and Vibration Research
%   Leibniz University Hannover
%   03.11.2020 - Tido Kubatschek
%   21.02.2021 - Alwin Förster
%
function [Do,Info,Path,Jacobian,breakFunOut,Opt,Counter,ds] = exitLoop(Do, Info, Initial, Is, Path, Opt, Counter, Bifurcation, ds, funSolution, Jacobian, breakFunOut)
    %% eval. break function:
    %
    try
        if Is.valid
            [bfun,breakFunOut] = Opt.breakFunction(funSolution,Jacobian.solver,Path.varAll,Path.lAll,breakFunOut);
        else
            bfun = false;
        end
    catch
        aux.printLine(Opt,'--> Unable to evaluate user defined break function.\n');
        bfun = false;
    end
    %
    %% exit because time limit reached:
    %
    if toc(Info.t0)>=Opt.timeLimit
        Info.exitflag = -5;
        Do.continuation = false;
        Info.exitMsg = sprintf('--> continuation stoped: Time limit was reached.');
    end
    %
    %% exit without complete results:
    %
    if Counter.catch>=3
        Info.exitflag = -4;
        Do.continuation = false;
        Info.exitMsg = sprintf('--> continuation stoped: Too many erros in function evaluation.');
    end
    %
    %% exit with maxRemoveCounter reached:
    %
    if Counter.remove>Opt.maxRemoveCounter
        Do.continuation = false;
        Info.exitflag = -3;
        Info.exitMsg = '--> continuation stoped: maxRemoveCounter has been reached.';
    end
    %
    %% exit without initial solution:
    %
    % exitflag = -2
    %
    %% exit without complete results:
    %
    if Counter.error>=Opt.maxErrorCounter
        Info.exitflag = -1;
        Do.continuation = false;
        Info.exitMsg = sprintf('--> continuation stoped: No valid result could be found for the last %d attempts.',Opt.maxErrorCounter);
    end
    %
    %% exit with l<lStart:
    %
    if sign(Info.lEnd-Info.lStart)*(Path.lAll(end)-Info.lStart)<0
        Do.continuation = false;
        Info.exitflag = 0;
        Info.exitMsg = '--> continuation stoped: l<lStart';
    end
    %
    %% exit with success:
    %
    if sign(Info.lEnd-Info.lStart)*(Path.lAll(end)-Info.lEnd)>=0 || sign(Opt.lTarget-Info.lStart)*(Path.lAll(end)-Opt.lTarget)>=0 || (Do.convergeToxTarget && norm(Opt.xTarget-[Path.varAll(:,end);Path.lAll(end)])<10^-8)
        Do.continuation = false;
        Info.exitflag = 1;
        if sign(Info.lEnd-Info.lStart)*(Path.lAll(end)-Info.lEnd)>=0
            Info.exitMsg = '--> continuation completed: lEnd reached';
        elseif sign(Opt.lTarget-Info.lStart)*(Path.lAll(end)-Opt.lTarget)>=0
            Info.exitMsg = '--> continuation completed: lTarget reached';
        elseif (Do.convergeToxTarget && norm(Opt.xTarget-[Path.varAll(:,end);Path.lAll(end)])<10^-8)
            Info.exitMsg = '--> continuation completed: xTarget reached';
        else
            Info.exitMsg = '--> continuation completed: unknown reason';
        end
    end
    %
    %% exit with nStepMax reached:
    %
    if Counter.loop>=Opt.nStepMax
        Do.continuation = false;
        Info.exitflag = 2;
        Info.exitMsg = '--> continuation completed: nStepMax reached';
    end
    %
    %% exit with bifurcation:
    %
    if Bifurcation.flag>0 && Opt.stopOnBifurcation
        Do.continuation = false;
        Info.exitflag = 3;
        Path.varAll = Path.varAll(:,1:Bifurcation.bif(1,end));
        Path.lAll = Path.lAll(1:Bifurcation.bif(1,end));
        Path.sAll = Path.sAll(1:Bifurcation.bif(1,end));
        if Opt.jacobianOut.full
            Jacobian.all(:,:,1:Bifurcation.bif(1,end))
        end
        Info.exitMsg = '--> continuation completed: bifurcation reached';
    end
    %
    %% exit on closed curve:
    %
    if Opt.closedCurveDetection
        [isClosed, Opt, Counter] = aux.closedCurve(Opt,Path,ds,Counter);
        if isClosed
            Do.continuation = false;
            Info.exitflag = 4;
            Info.exitMsg = '--> continuation completed: closed curve detected';
        end
    end
    %
    %% exit due to break function:
    %
    if bfun
        Do.continuation = false;
        Info.exitflag = 5;
        Info.exitMsg = '--> continuation completed: user-defined break function';
    end
    %
    %% exit due to user input:
    %
    if Do.stopManually
        Do.continuation = false;
        Info.exitflag = 6;
        Info.exitMsg = '--> continuation stoped by user';
    end
    %
    %% exit with bifurcation:
    %
    if Bifurcation.flag>0 && Opt.stopOnCrossing && Bifurcation.bif(2,end)==0
        Do.continuation = false;
        Info.exitflag = 7;
        Path.varAll = Path.varAll(:,1:Bifurcation.bif(1,end));
        Path.lAll = Path.lAll(1:Bifurcation.bif(1,end));
        Path.sAll = Path.sAll(1:Bifurcation.bif(1,end));
        if Opt.jacobianOut.full
            Jacobian.all(:,:,1:Bifurcation.bif(1,end))
        end
        Info.exitMsg = '--> continuation completed: bifurcation reached';
    end
    %
    %% check bidirectional runs:
    %
    if ~Do.continuation && Opt.bidirectional && (Info.biDirRuns==0)
        Info.biDirRuns = Info.biDirRuns+1;
        Do.continuation = true;
        % turn path
        Path.varAll = Path.varAll(:,end:-1:1);
        Path.lAll = Path.lAll(end:-1:1);
        Path.sAll = Path.sAll(end)-Path.sAll(end:-1:1);
        if Opt.jacobianOut.full
            Jacobian.all = Jacobian.all(:,:,end:-1:1);
        end
        if Opt.lTarget==Info.lEnd
            Opt.lTarget = Info.lStart;
            Opt.direction = -Opt.direction;
            lStartTemp = Info.lStart;
            Info.lStart = Info.lEnd;
            Info.lEnd = lStartTemp;
        end
        % reset values
        Path.speedOfContinuation = [];
        Path.xPredictor = [];
        Path.bifTestValue = [];
        Path.pathInfoValue = [];
        Jacobian.last = [];
        Jacobian.previous = [];
        Jacobian.signDet = sign(det(Jacobian.initial));
        Jacobian.solver = [];
        % set step size
        if numel(Path.sAll)>1
            ds = Path.sAll(end)-Path.sAll(end-1);
        else
            ds = Info.ds0;
        end
    elseif ~Do.continuation && Opt.bidirectional && (Info.biDirRuns>0)
        if Initial.lStart~=Info.lStart
            % turn path
            Path.varAll = Path.varAll(:,end:-1:1);
            Path.lAll = Path.lAll(end:-1:1);
            Path.sAll = Path.sAll(end)-Path.sAll(end:-1:1);
            if Opt.jacobianOut.full
                Jacobian.all = Jacobian.all(:,:,end:-1:1);
            end
            % reset values
            Path.speedOfContinuation = [];
            Path.xPredictor = [];
            Path.bifTestValue = [];
            Path.pathInfoValue = [];
            Jacobian.last = [];
            Jacobian.previous = [];
            Jacobian.signDet = sign(det(Jacobian.initial));
            Jacobian.solver = [];
        end
    end
    %
    %% reset
    %
    Do.convergeToxTarget = false;
    %
end