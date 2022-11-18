%% path continuation - continuation
%
%   Institute of Dynamics and Vibration Research
%   Leibniz University Hannover
%   08.05.2020 - Alwin F�rster
%
%   [varAll,lAll,exitflag,Bifurcation] = continuation(fun,var0,lStart,lEnd,ds0,varargin)
%
%   fun = fun(var,l) != 0
%   lStart <= l <= lEnd
%   ds0: initial stepsize
%
%% This file is part of continuation.
% 
% If you use continuation, please refer to:
%   A. F�rster, foerster@ids.uni-hannover.de
% 
% COPYRIGHT AND LICENSING: 
% Continuation Copyright (C) 2022 Alwin F�rster
%                                 (foerster@ids.uni-hannover.de)
%                                 Leibnitz University Hannover
% This program comes with NO WARRANTY. 
% Continuation is free software, you can redistribute and/or modify it
% under the GNU General Public License as published by the Free Software
% Foundation, either version 3 of the License, or (at your option) any
% later version. For details on license and warranty, see
% http://www.gnu.org/licenses or gpl-3.0.txt.
%
%%
function [varAll,lAll,exitflag,Bifurcation,sAll,jacobianOut,breakFunOut,InfoOut] = ...
    continuation(fun,var0,lStart,lEnd,ds0,varargin)
    %% initialize
    %
    warning on;
    [Opt,ds0,OptIsSet,func] = continuation.input(varargin,fun,var0,lStart,lEnd,ds0);
    [Opt,ds0,StepsizeOptions] = stepSize.initialize(Opt,var0,lStart,lEnd,ds0);
    if StepsizeOptions.rateOfContraction
        global solverStepsizes;
    end
    [Bifurcation,Counter,Do,Info,InfoOut,Initial,Is,Jacobian,Path,Plot,Plus,Remove,Solver,StepsizeInformation,Temp] = aux.initializeStructs(var0,lStart,lEnd,ds0,Opt,StepsizeOptions.rateOfContraction);
    clear('var0','lStart','lEnd','ds0');
    resCorr = continuation.corrector(fun,Opt);
    ds = Info.ds0;
    aux.printLine(Opt,'Starting path continuation...\n');
    if Opt.dpaGammaVar
        paraName = 'g';
    else
        paraName = 'l';
    end
    tDisplay = tic;
    %
    %% find initial solution
    %
    residualInitial = @(v) aux.residualFixedValue(func,v,Opt.l0,Opt);
    [Path.varAll,funInitial,initialExitflag,Solver.output,Jacobian.initial] = Solver.main(residualInitial,Info.var0,Opt.dscale0(1:end-1));
    Jacobian.solver = Jacobian.initial;
    breakFunOut = [];
    event.eventObj = [];
    if initialExitflag>=0
        Path.lAll = Opt.l0;
        Path.sAll = 0;
        Path.biftestValue = Opt.bifAdditionalTestfunction(func,[Path.varAll;Path.lAll],Jacobian,Path,Info);
        Do.continuation = true;
        Do.loop = true;
        dpaPoints = [];
        Jacobian.solver = [Jacobian.solver,aux.numericJacobian(@(x) func(x(1:Info.nv),x(Info.nv+1)),[Path.varAll;Opt.l0],'centralValue',funInitial,'derivativeDimensions',Info.nv+1,'diffquot',Opt.diffquot)];
        Jacobian.previous = Jacobian.solver;
        Jacobian.last = Jacobian.solver;
        [~,breakFunOut] = Opt.breakFunction(funInitial,Jacobian.solver,Path.varAll,Path.lAll,breakFunOut);
        aux.printLine(Opt,'Initial solution at %s = %.2e\n',paraName,Opt.l0);
        if aux.ison(Opt.bifurcation)
            Jacobian.signDetRed = sign(det(Jacobian.initial));
        end
        if aux.ison(Opt.plot)
            [Plot, Opt] = plot.livePlot(Opt, Info, Path, Info.ds0, Info.ds0, Solver.output.iterations(end), Counter);
        end
        if StepsizeOptions.rateOfContraction
            if size(solverStepsizes, 1) < 3
                Solver.output.rateOfContraction = Opt.optimalContractionRate;
            else
                Solver.output.rateOfContraction = solverStepsizes(3,2)/solverStepsizes(2,2);
            end
        end
        if StepsizeOptions.speedOfContinuation
             Path.speedOfContinuation = Opt.speedOfContinuation;
        end
        %
        if StepsizeOptions.predictor
            Path.xPredictor = [Info.var0;Info.lStart];
        end
    else
        Path.varAll = [];
        Path.lAll = [];
        Path.sAll = [];
        Info.exitflag = -2;
        Do.continuation = false;
        aux.printLine(Opt,'No initial solution found.\n');
    end
    %
    %% continuation
    %
    while Do.continuation
        %% initialize loop
        %
        Counter.loop = Counter.loop+1;
        Is.currentJacobian = false;
        Counter.catchOld = Counter.catch;
        Temp = stepSize.updateTemp(Path,Solver,StepsizeOptions,Temp);
        %
        %% residual
        %
        if Do.deflate
            try
                residual = @(x) aux.deflation(residual,xDeflation,x,Opt);
            catch
                aux.printLine(Opt,'---> delation: catch!\n');
                Counter.catch = Counter.catch + 1;
            end
        elseif Do.suspend
            residual = @(v,lFix) aux.residualSuspendContinuation(func,v,lFix,Opt);
        else
            residual = @(x) aux.mergeResiduals(func,resCorr,x,[Path.varAll;Path.lAll],ds,Jacobian.last,Opt);
        end
        if numel(Path.lAll)==1 && ~isempty(Opt.initialDeflationPoints)
            for ii=1:numel(Opt.initialDeflationPoints(1,:))
                residual = @(x) aux.deflation(residual,Opt.initialDeflationPoints(:,ii),x,Opt);
            end
        end
        %
        %% predictor
        %
        if StepsizeOptions.speedOfContinuation
            timeNeeded = tic;
        end
        %
        try
            if Do.homotopy
                %% Homotopy
                xPredictor = homotopy.hContinuation(residual,[Path.varAll(:,end);Path.lAll(end)],Opt);
                varPredictor = xPredictor(1:end-1);
                lPredictor = xPredictor(end);
            else
                %% calc. predictor
                [varPredictor,lPredictor,funPredictor,sPredictor,ds] = continuation.predictor(Path,ds,Jacobian.last,func,resCorr,Solver,Opt);
                xPredictor = [varPredictor;lPredictor];
            end
        catch
            [varPredictor,lPredictor,funPredictor,sPredictor,ds] = continuation.predictor(Path,ds,Jacobian.last,func,resCorr,Solver,Opt);
            xPredictor = [varPredictor;lPredictor];
            aux.printLine(Opt,'---> predictor: catch!\n');
            Counter.catch = Counter.catch + 1;
        end
        %
        %% solve
        %
        try
            dscale = aux.getDscale(Opt,Path);
            if sign(Path.lAll(end)-Opt.lTarget)*sign(lPredictor-Opt.lTarget)<=0
                %% try to converge to target
                residualTarget = @(v) aux.residualFixedValue(func,v,Opt.lTarget,Opt);
                varPredictorCtt = (varPredictor - Path.varAll(:,end))*(abs(Opt.lTarget-Path.lAll(end))/abs(lPredictor-Path.lAll(end)))+Path.varAll(:,end);
                [varSolution,funSolution,Solver.exitflag,Solver.output,Jacobian.solver] = Solver.main(residualTarget,varPredictorCtt,dscale(1:end-1));
                xSolution = [varSolution;Opt.lTarget];
                Do.convergeToTarget = true;
            else
                %% regular solver
                %            
                if Opt.solver.fsolve && Opt.solverForce1it
                    [xSolution,funSolution,Solver.exitflag,Solver.output,Jacobian.solver] = Solver.main(residual,xPredictor,dscale);
                    if Solver.output.iterations(end) < 1
                        % perturbate initial solution by tolerance of
                        % solver
                        pert = Opt.solverTol * ones(numel(xPredictor),1) / numel(xPredictor);
                        [xSolution,funSolution,Solver.exitflag,Solver.output,Jacobian.solver] = Solver.main(residual,xPredictor + pert,dscale);
                    end
                else
                    if Do.suspend
                        [vSolution,funSolution,Solver.exitflag,Solver.output,Jacobian.solver] = Solver.main(@(v) residual(v,xPredictor(end)),xPredictor(1:(end-1)),dscale(1:(end-1)));
                        xSolution = [vSolution;xPredictor(end)];
                    else
                        [xSolution,funSolution,Solver.exitflag,Solver.output,Jacobian.solver] = Solver.main(residual,xPredictor,dscale);
                    end
                end
                Do.convergeToTarget = false;
            end
            Is.currentJacobian = true;
        catch
            xSolution = NaN(size(xPredictor));
            funSolution = inf(size(xPredictor));
            Solver.exitflag = -2;
            Solver.output = Solver.defaultOutput;
            Do.convergeToTarget = false;
            aux.printLine(Opt,'---> solve: catch!\n');
            Counter.catch = Counter.catch + 1;
        end
        %
        %% calc stepsize information
        %
        % measure speed
        %
        if StepsizeOptions.speedOfContinuation
            timeNeeded = toc(timeNeeded);
            StepsizeInformation.speedOfContinuation = ds/timeNeeded; 
        end
        %
        % measure rate of contraction
        if StepsizeOptions.rateOfContraction
            if size(solverStepsizes, 1) < 3
                StepsizeInformation.rateOfContraction = Opt.optimalContractionRate;
            else
                StepsizeInformation.rateOfContraction = solverStepsizes(3,2)/solverStepsizes(2,2);
            end
        end
        %
        if StepsizeOptions.iterations
            if ~isempty(Temp.iterations)
                Temp.iterations = [Temp.iterations, Solver.output.iterations];
            else
                Temp.iterations = Solver.output.iterations;
            end
        end
        %
        %% adaptive corrector
        %
        [Do,Opt,corrInfo] = corrector.adapt(Do,Opt,Path,Solver,func,xPredictor,dscale,Jacobian.last,ds);
        %
        %% check result
        %
        % check result:
        [invPoiStr,Counter,Do,Is,Opt] = ...
            aux.validateResult(xSolution,Plus,funSolution,Path,ds,Solver,Jacobian,funPredictor,sPredictor,Do,Bifurcation,Info,Is,Counter,Plot,Opt);
        % confirm result:
        [xDeflation,Bifurcation,Counter,Do,Info,Initial,Is,Jacobian,Path,Plus,Remove,Solver,StepsizeInformation,StepsizeOptions,Temp,Opt] = ...
            aux.confirmResult(func,xSolution,xPredictor,Bifurcation,Counter,Do,Info,Initial,Is,Jacobian,Path,Plus,Remove,Solver,StepsizeInformation,StepsizeOptions,Temp,Opt,OptIsSet);
        %
        %% Bifurcations
        %
        Bifurcation.flag = 0;
        if aux.ison(Opt.bifurcation) && Is.valid && ~Do.homotopy && numel(Path.lAll)>2
            if ~Is.currentJacobian
                %% get jacobian if not current
                Jacobian.solver = aux.getJacobian(func,Path.varAll(:,end),Path.lAll(end),Opt);
            end
            [Bifurcation,Jacobian,Path] = bifurcation.check(func,Jacobian,Path,Bifurcation,Info,resCorr,Solver,Opt,OptIsSet);
        elseif aux.ison(Opt.bifurcation) && Is.valid && numel(Path.lAll)<=2
            if ~Is.currentJacobian
                %% get jacobian if not current
                Jacobian.solver = aux.getJacobian(func,Path.varAll(:,end),Path.lAll(end),Opt);
            end
            Jacobian.signDetRed = sign(det(Jacobian.solver(1:Info.nv,1:Info.nv)));
        end
        %
        %% DPA points
        %
        if Opt.dpa && OptIsSet.dpa && ~Opt.dpaGammaVar
            [Path,dpaPoints] = dpa.checkResidual(fun,dpaPoints,Opt,Path,Solver);
        end
        %
        %% step size control
        %
        % save latest stepsize:
        dsim1 = ds;
        % save step size data:
        [Solver,Path] = aux.updateStepsizeData(StepsizeOptions,Temp,Solver,Path);
        % adjust stepsize:
        [ds,Counter,event,Opt] = stepSize.control(ds,Counter,Solver,Do,Plus,Path,Jacobian,Opt,Info,event,Initial);
        %
        %% end loop
        %
        if Is.valid
            if ~isempty(Solver.output.iterations)
                aux.printLine(Opt,'-----> continued at %s = %.4e\t|\tnew step size: ds = %.2e\t|\tloop counter = %d\t|\tstep = %d\t|\titerations = %d/%d\n',paraName,Path.lAll(end),ds,Counter.loop,Counter.step,Solver.output.iterations(end),Opt.nIterOpt);
            else
                aux.printLine(Opt,'-----> continued at %s = %.4e\t|\tnew step size: ds = %.2e\t|\tloop counter = %d\t|\tstep = %d\t|\titerations = %d/%d\n',paraName,Path.lAll(end),ds,Counter.loop,Counter.step,[],Opt.nIterOpt);
            end
        else
            if ~isempty(Solver.output.iterations)
                aux.printLine(Opt,'-----> invalid point %s |\tnew step size: ds = %.2e\t|\tloop counter = %d\t|\tstep = %d\t|\titerations = %d/%d\n',invPoiStr,ds,Counter.loop,Counter.step,Solver.output.iterations(end),Opt.nIterOpt);
            else
                aux.printLine(Opt,'-----> invalid point %s |\tnew step size: ds = %.2e\t|\tloop counter = %d\t|\tstep = %d\t|\titerations = %d/%d\n',invPoiStr,ds,Counter.loop,Counter.step,[],Opt.nIterOpt);
            end
        end
        [Do,Info,Path,breakFunOut,Opt,Counter] = aux.exitLoop(Do,Info,Is,Path,Opt,Counter,Bifurcation,ds,funSolution,Jacobian,breakFunOut);
        if Do.changeCorrector
            Opt = aux.seton(Opt,'corrector',corrInfo);
            resCorr = continuation.corrector(fun,Opt);
            Do.changeCorrector = false;
        end
        Opt = aux.updateOpt(Opt,OptIsSet,Info);
        %
        %% live plot
        %
        if aux.ison(Opt.plot) && Is.valid
            try
                [Plot, Opt] = plot.livePlot(Opt, Info, Path, ds, dsim1, Solver.output.iterations(end), Counter, funPredictor, sPredictor, Plot, Bifurcation, dpaPoints);
            catch
                aux.printLine(Opt,'--> The plot update has failed.\n');
            end
        end
        %
        %% catch counter
        %
        if Counter.catchOld == Counter.catch
            Counter.catch = 0;
        end
        %
    end
    %
    %% bifurcation tracing
    %
    if Opt.bifurcation.trace
        try
            delete(Plot.plCurr);
            [Path,Bifurcation] = bifurcation.trace(Opt,OptIsSet,Path,Bifurcation,Solver,Info,func,resCorr);
            Jacobian.last = [];
        catch
            aux.printLine(Opt,'--> Failed to trace bifurcations.\n');
        end
    elseif Opt.bifurcation.parameterTrace
        try
            delete(Plot.plCurr);
            Path = bifurcation.parameterTrace(Opt,Path,Bifurcation,Info,fun);
            Jacobian.last = [];
        catch
            aux.printLine(Opt,'--> Failed to trace bifurcation parameter.\n');
        end
    end
    %
    %% DPA
    %
    if Opt.dpa && OptIsSet.dpa && ~Opt.dpaGammaVar && ~Opt.bifurcation.parameterTrace
        Path = dpa.trace(fun,dpaPoints,Info,Opt,Path);
    end
    %
    %% live plot finalization
    %
    if aux.ison(Opt.plot) && initialExitflag>0
        try
            BifurcationLastPlot = Bifurcation;
            BifurcationLastPlot.flag = -1;
            plot.livePlot(Opt, Info, Path, ds, dsim1, Solver.output.iterations(end), Counter, funPredictor, sPredictor, Plot, BifurcationLastPlot,dpaPoints);
            if isfield(Plot,'plCurr')
                delete(Plot.plCurr);
            end
        catch
            aux.printLine(Opt,'--> The plot update has failed.\n');
        end
        drawnow;
    end
    %
    %% final disp
    %
    aux.printLine(Opt,[Info.exitMsg,'\n']);
    aux.printLine(Opt,'--> time elapsed: %.3f s\n',toc(tDisplay));
    %
    %% output
    %
    InfoOut.numberOfSteps = Counter.step;
    InfoOut.numberOfInvalidPoints = Counter.loop - Counter.step;
    jacobianOut = Jacobian.last;
    exitflag = Info.exitflag;
    varAll = Path.varAll;
    lAll = Path.lAll;
    sAll = Path.sAll;
	%
end