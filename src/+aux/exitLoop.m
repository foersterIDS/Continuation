%% path continuation - aux.exitLoop
%
%   Institute of Dynamics and Vibration Research
%   Leibniz University Hannover
%   03.11.2020 - Tido Kubatschek
%   21.02.2021 - Alwin Förster
%
function [Do, Info, Path, breakFunOut, Opt,Counter] = exitLoop(Do, Info, Is, Path, Opt, Counter, Bifurcation, ds, funSolution, Jacobian, breakFunOut)
    %% eval. break function:
    %
    try
        if Is.valid
            [bfun,breakFunOut] = Opt.breakFunction(funSolution,Jacobian.solver,Path.varAll(:,end),Path.lAll(end),breakFunOut);
        else
            bfun = false;
        end
    catch
        aux.printLine(Opt,'--> Unable to evaluate user defined break function.\n');
        bfun = false;
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
    if sign(Info.lEnd-Info.lStart)*(Path.lAll(end)-Info.lEnd)>=0 || sign(Opt.lTarget-Info.lStart)*(Path.lAll(end)-Opt.lTarget)>=0
        Do.continuation = false;
        Info.exitflag = 1;
        if sign(Info.lEnd-Info.lStart)*(Path.lAll(end)-Info.lEnd)>=0
            Info.exitMsg = '--> continuation completed: lEnd reached';
        else
            Info.exitMsg = '--> continuation completed: lTarget reached';
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
        Info.exitMsg = '--> continuation completed: bifurcation reached';
    end
    %
end

