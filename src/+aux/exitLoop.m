%% path continuation - aux.exitLoop
%
%   Institute of Dynamics and Vibration Research
%   Leibniz University Hannover
%   03.11.2020 - Tido Kubatschek
%   21.02.2021 - Alwin Förster
%
function [breakFunOut,ds] = exitLoop(oih, ds, funSolution, breakFunOut)
    %% eval. break function:
    %
    try
        if oih.is.valid && oih.optIsSet.breakFunction
            oih.path.outputFormat = 'full';
            [bfun,breakFunOut] = oih.opt.breakFunction(funSolution,oih.solver.jacobian,oih.path.varAll,oih.path.lAll,breakFunOut);
            oih.path.resetOutput();
        else
            bfun = false;
        end
    catch
        aux.printLine(oih,'--> Unable to evaluate user defined break function.\n');
        bfun = false;
    end
    %
    %% exit because time limit reached:
    %
    if toc(oih.info.t0)>=oih.opt.timeLimit
        oih.info.exitflag = -5;
        oih.do.continuation = false;
        oih.info.exitMsg = sprintf('--> continuation stoped: Time limit was reached.');
    end
    %
    %% exit without complete results:
    %
    if oih.counter.catch>=3
        oih.info.exitflag = -4;
        oih.do.continuation = false;
        oih.info.exitMsg = sprintf('--> continuation stoped: Too many erros in function evaluation.');
    end
    %
    %% exit with maxRemoveCounter reached:
    %
    if oih.counter.remove>oih.opt.maxRemoveCounter
        oih.do.continuation = false;
        oih.info.exitflag = -3;
        oih.info.exitMsg = '--> continuation stoped: maxRemoveCounter has been reached.';
    end
    %
    %% exit without initial solution:
    %
    % exitflag = -2
    %
    %% exit without complete results:
    %
    if oih.counter.error>=oih.opt.maxErrorCounter
        oih.info.exitflag = -1;
        oih.do.continuation = false;
        oih.info.exitMsg = sprintf('--> continuation stoped: No valid result could be found for the last %d attempts.',oih.opt.maxErrorCounter);
    end
    %
    %% exit with l<lStart:
    %
    if sign(oih.info.lEnd-oih.info.lStart)*(oih.path.lAll(end)-oih.info.lStart)<=0
        oih.do.continuation = false;
        oih.info.exitflag = 0;
        oih.info.exitMsg = '--> continuation stoped: l<=lStart';
    end
    %
    %% exit with success:
    %
    if sign(oih.info.lEnd-oih.info.lStart)*(oih.path.lAll(end)-oih.info.lEnd)>=0 || sign(oih.opt.lTarget-oih.info.lStart)*(oih.path.lAll(end)-oih.opt.lTarget)>=0 || (oih.do.convergeToxTarget && norm(oih.opt.xTarget-[oih.path.varAll(:,end);oih.path.lAll(end)])<10^-8)
        oih.do.continuation = false;
        oih.info.exitflag = 1;
        if sign(oih.info.lEnd-oih.info.lStart)*(oih.path.lAll(end)-oih.info.lEnd)>=0
            oih.info.exitMsg = '--> continuation completed: lEnd reached';
        elseif sign(oih.opt.lTarget-oih.info.lStart)*(oih.path.lAll(end)-oih.opt.lTarget)>=0
            oih.info.exitMsg = '--> continuation completed: lTarget reached';
        elseif (oih.do.convergeToxTarget && norm(oih.opt.xTarget-[oih.path.varAll(:,end);oih.path.lAll(end)])<10^-8)
            oih.info.exitMsg = '--> continuation completed: xTarget reached';
        else
            oih.info.exitMsg = '--> continuation completed: unknown reason';
        end
    end
    %
    %% exit with nStepMax reached:
    %
    if oih.counter.loop>=oih.opt.nStepMax
        oih.do.continuation = false;
        oih.info.exitflag = 2;
        oih.info.exitMsg = '--> continuation completed: nStepMax reached';
    end
    %
    %% exit with bifurcation:
    %
    if oih.bifurcation.flag>0 && oih.opt.stopOnBifurcation
        oih.do.continuation = false;
        oih.info.exitflag = 3;
        oih.path.remove((oih.bifurcation.bif(1,end)+1):oih.path.nAll);
        oih.info.exitMsg = '--> continuation completed: bifurcation reached';
    end
    %
    %% exit on closed curve:
    %
    if oih.opt.closedCurveDetection
        [isClosed] = aux.closedCurve(oih,ds);
        if isClosed
            oih.do.continuation = false;
            oih.info.exitflag = 4;
            oih.info.exitMsg = '--> continuation completed: closed curve detected';
        end
    end
    %
    %% exit due to break function:
    %
    if bfun
        oih.do.continuation = false;
        oih.info.exitflag = 5;
        oih.info.exitMsg = '--> continuation completed: user-defined break function';
    end
    %
    %% exit due to user input:
    %
    if oih.do.stopManually
        oih.do.continuation = false;
        oih.info.exitflag = 6;
        oih.info.exitMsg = '--> continuation stoped by user';
    end
    %
    %% exit with bifurcation:
    %
    if oih.bifurcation.flag>0 && oih.opt.stopOnCrossing && oih.bifurcation.bif(2,end)==0
        oih.do.continuation = false;
        oih.info.exitflag = 7;
        oih.path.remove((oih.bifurcation.bif(1,end)+1):oih.path.nAll);
        oih.info.exitMsg = '--> continuation completed: bifurcation reached';
    end
    %
    %% check bidirectional runs:
    %
    if ~oih.do.continuation && oih.opt.bidirectional && (oih.info.biDirRuns==0)
        oih.info.biDirRuns = oih.info.biDirRuns+1;
        oih.do.continuation = true;
        % turn path
        oih.path.turn();
        if oih.opt.lTarget==oih.info.lEnd
            oih.opt.lTarget = oih.info.lStart;
            oih.opt.direction = -oih.opt.direction;
            lStartTemp = oih.info.lStart;
            oih.info.lStart = oih.info.lEnd;
            oih.info.lEnd = lStartTemp;
        end
        % reset values
        oih.solver.jacobian = [];
        % set step size
        if numel(oih.path.sAll)>1
            ds = oih.path.sAll(end)-oih.path.sAll(end-1);
        else
            ds = oih.info.ds0;
        end
    elseif ~oih.do.continuation && oih.opt.bidirectional && (oih.info.biDirRuns>0)
        if oih.initial.lStart~=oih.info.lStart
            % turn path
            oih.path.turn();
            % reset values
            oih.solver.jacobian = [];
        end
    end
    %
    %% reset
    %
    oih.do.convergeToxTarget = false;
    %
end