%% path continuation - continuation
%
%   Institute of Dynamics and Vibration Research
%   Leibniz University Hannover
%   08.05.2020 - Alwin Förster
%
%   [varAll,lAll,exitflag,Bifurcation] = continuation(fun,var0,lStart,lEnd,ds0,NameValueArgs)
%
%   fun = fun(var,l) != 0
%   lStart <= l <= lEnd
%   ds0: initial stepsize
%
%% This file is part of continuation.
% 
% If you use continuation, please refer to:
%   A. Förster, foerster@ids.uni-hannover.de
% 
% COPYRIGHT AND LICENSING: 
% Continuation Copyright (C) 2023 Alwin Förster
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
    continuation(fun,var0,lStart,lEnd,ds0,NameValueArgs)
    %% arguments
    %
    arguments
        fun (1,1) function_handle
        var0 (:,1) double
        lStart (1,1) double
        lEnd (1,1) double
        ds0 (1,1) double {mustBeGreaterThan(ds0,0)}
        NameValueArgs.Opt (1,1) struct
        NameValueArgs.adaptCorrector (1,:) char {mustBeMember(NameValueArgs.adaptCorrector,{'basic','solve'})}
        NameValueArgs.alphaReverse (1,1) double {mustBeGreaterThan(NameValueArgs.alphaReverse,0)}
        NameValueArgs.approveManually {validation.scalarLogical} % #scalar#logical#ison:plot#ison:display|#scalar#double#ison:plot#ison:display
        NameValueArgs.bidirectional {validation.scalarLogical}
        NameValueArgs.bifurcation (1,:) char {mustBeMember(NameValueArgs.bifurcation,{'on','off','mark','determine','trace','parameterTrace'})}
        NameValueArgs.bifRandDir {validation.scalarLogical}
        NameValueArgs.bifAdditionalTestfunction (1,1) function_handle
        NameValueArgs.bifResidual (1,:) char {mustBeMember(NameValueArgs.bifResidual,{'determinant','luFactorization'})}
        NameValueArgs.breakFunction (1,1) function_handle
        NameValueArgs.checkResidual {validation.scalarLogical}
        NameValueArgs.checkStepSizeOptions {validation.scalarLogical}
        NameValueArgs.closedCurveDetection {validation.scalarLogical}
        NameValueArgs.corrector (1,:) char {mustBeMember(NameValueArgs.corrector,{'sphere','orthogonal','orthogonal2','ellipsoid','ellipsoid2','unique','paraboloid'})}
        NameValueArgs.correctorOrthogonalMethod (1,:) char {mustBeMember(NameValueArgs.correctorOrthogonalMethod,{'secant','tangent'})}
        NameValueArgs.correctPredictor {validation.scalarLogical}
        NameValueArgs.deflation {validation.scalarLogical}
        NameValueArgs.deflationErrorCounter (1,1) double {mustBeGreaterThan(NameValueArgs.deflationErrorCounter,0)}
        NameValueArgs.diffquot (1,:) char {mustBeMember(NameValueArgs.diffquot,{'forward','central'})}
        NameValueArgs.direction (:,1) double % #scalar#pmone|#array#double#norm:1#size:[numel(var0)+1,1]
        NameValueArgs.display {validation.scalarLogical}
        NameValueArgs.dpa {validation.scalarLogical} % #scalar#logical#false|#scalar#logical#true#isoff:plot|#scalar#logical#true#ison:plot.dpa
        NameValueArgs.dpaGammaVar {validation.scalarLogical} % #scalar#logical#false|#scalar#logical#true#ison:dpa#isoff:plot|#scalar#logical#true#ison:dpa#ison:plot.dpa
        NameValueArgs.dpaResidual (1,1) function_handle % #nargin:3
        NameValueArgs.dsMax (:,1) double {mustBeGreaterThan(NameValueArgs.dsMax,0)} % #scalar#double#positive#nonzero|#array#positive#nonzero#double#size:[numel(var0)+1,1]#ison:enforceDsMax
        NameValueArgs.dsMin (1,1) double {mustBeGreaterThan(NameValueArgs.dsMin,0)} % #scalar#double#positive#smaller:norm(Opt.dsMax)
        NameValueArgs.dsTol (1,2) double {mustBePositive} % #array#positive#double#increasing#size:[1,2]
        NameValueArgs.dscale0 (:,1) double {mustBeGreaterThan(NameValueArgs.dscale0,0)} % #array#positive#nonzero#double#size:[numel(var0)+1,1]
        NameValueArgs.dscaleMin (:,1) double {mustBeGreaterThan(NameValueArgs.dscaleMin,0)} % #scalar#positive#nonzero#double|#array#positive#nonzero#double#size:[numel(var0)+1,1]
        NameValueArgs.enforceDsMax {validation.scalarLogical} % #ison:predictorSolver
        NameValueArgs.eventUserInput (1,:) cell
        NameValueArgs.g0 (1,1) double
        NameValueArgs.gTarget (1,1) double
        NameValueArgs.homotopy (1,:) char {mustBeMember(NameValueArgs.homotopy,{'on','off','f2','fix','fixnt','newton','squared'})}
        NameValueArgs.homotopyErrorCounter (1,1) double {mustBeInteger,mustBeGreaterThan(NameValueArgs.homotopyErrorCounter,0)} % #neq:Opt.deflationErrorCounter
        NameValueArgs.includeReverse {validation.scalarLogical}
        NameValueArgs.initialDeflationPoints (1,1) double
        NameValueArgs.jacobian {validation.scalarLogical}
        NameValueArgs.l0 (1,1) double
        NameValueArgs.lTarget (1,1) double
        NameValueArgs.livePlotFig (1,1) double % #scalar#isnan|#scalar#integer#positive#nonzero
        NameValueArgs.maxClosedCounter (1,1) double {mustBeGreaterThan(NameValueArgs.maxClosedCounter,0)}
        NameValueArgs.maxErrorCounter (1,1) double {mustBeGreaterThan(NameValueArgs.maxErrorCounter,0)}
        NameValueArgs.maxRemoveCounter (1,1) double {mustBeGreaterThan(NameValueArgs.maxRemoveCounter,0)}
        NameValueArgs.maxStepSizeChange (1,1) double {mustBeGreaterThan(NameValueArgs.maxStepSizeChange,1)}
        NameValueArgs.nBifSearch (1,1) double {mustBeGreaterThan(NameValueArgs.nBifSearch,0)}
        NameValueArgs.nIterOpt (1,1) double {mustBeGreaterThan(NameValueArgs.nIterOpt,0)}
        NameValueArgs.nStepMax (1,1) double {mustBeInteger,mustBeGreaterThan(NameValueArgs.nStepMax,0)}
        NameValueArgs.optimalContractionRate (1,1) double {mustBeGreaterThan(NameValueArgs.optimalContractionRate,0),mustBeSmallerThan(NameValueArgs.optimalContractionRate,1)}
        NameValueArgs.pathInfoFunction (1,1) function_handle
        NameValueArgs.plot (1,:) char {mustBeMember(NameValueArgs.plot,{'on','off','basic','detail','dpa','semilogx','semilogy','loglog','threeDim'})}
        NameValueArgs.plotPause {validation.scalarLogical} % #scalar#positive#nonzero#integer|#scalar#logical
        NameValueArgs.plotVarOfInterest (1,1) double {mustBeGreaterThan(NameValueArgs.plotVarOfInterest,0)} % #scalar#isnan|#scalar#integer#positive#nonzero#max:numel(var0)
        NameValueArgs.plotVarsIndex (1,:) double % #array#integer#positive#nonzero#unique#max:numel(var0)#ison:plot
        NameValueArgs.predictor (1,:) char {mustBeMember(NameValueArgs.predictor,{'polynomial','tangential'})}
        NameValueArgs.predictorDistance (1,1) double {mustBeGreaterThan(NameValueArgs.predictorDistance,0)}
        NameValueArgs.predictorPolynomialAdaptive {validation.scalarLogical}
        NameValueArgs.predictorPolynomialFit (1,1) double {mustBeInteger,mustBePositive}
        NameValueArgs.predictorPolynomialDegree (1,1) double {mustBeGreaterThan(NameValueArgs.predictorPolynomialDegree,0)}
        NameValueArgs.predictorSolver {validation.scalarLogical}
        NameValueArgs.removeErrorCounter (1,1) double {mustBeInteger} % #scalar#integer#positive#nonzero#larger:Opt.stepbackErrorCounter+1#neq:Opt.homotopyErrorCounter#neq:Opt.deflationErrorCounter#neq:Opt.suspendContinuationErrorCounter|#scalar#integer#equals:0
        NameValueArgs.reverse  (1,:) char {mustBeMember(NameValueArgs.reverse,{'angle','jacobian'})}
        NameValueArgs.scaling (1,:) char {mustBeMember(NameValueArgs.scaling,{'dynamicdscale','staticdscale'})}
        NameValueArgs.solver (1,:) char {mustBeMember(NameValueArgs.solver,{'fsolve','lsqnonlin','newton'})}
        NameValueArgs.solverForce1it {validation.scalarLogical}
        NameValueArgs.solverMaxIterations (1,1) double {mustBeGreaterThan(NameValueArgs.solverMaxIterations,0),mustBeInteger}
        NameValueArgs.solverTol (1,1) double {mustBeGreaterThan(NameValueArgs.solverTol,0)}
        NameValueArgs.speedOfContinuation (1,1) double {mustBeGreaterThan(NameValueArgs.speedOfContinuation,0)}
        NameValueArgs.stepbackErrorCounter (1,1) double {mustBePositive,mustBeInteger,mustBeGreaterThan(NameValueArgs.stepbackErrorCounter,0)} % #neq:Opt.homotopyErrorCounter#neq:Opt.deflationErrorCounter
        NameValueArgs.stepSizeAngle (1,1) double {mustBePositive}
        NameValueArgs.stepSizeControl (1,:) char {mustBeMember(NameValueArgs.stepSizeControl,{'angleChange','angleCustom','contraction','error','errorAlt','fayezioghani','fix','iterationsExponential','iterationsPolynomial','multiplicative','multiplicativeAlt','pidCustom','pidValli','szyszkowski','yoon'})}
        NameValueArgs.stepSizeErrorMax (1,1) double {mustBeGreaterThan(NameValueArgs.stepSizeErrorMax,0)}
        NameValueArgs.stepSizeErrorPd (1,:) double {mustBePositive}
        NameValueArgs.stepSizeEvent {validation.scalarLogical}
        NameValueArgs.stepSizeIterationsBeta (1,1) double {mustBeGreaterThan(NameValueArgs.stepSizeIterationsBeta,0),mustBeSmallerThan(NameValueArgs.stepSizeIterationsBeta,2)}
        NameValueArgs.stepSizeExponentialWeight (1,1) double {mustBeGreaterThan(NameValueArgs.stepSizeExponentialWeight,0)}
        NameValueArgs.stepSizePidParams (1,3) double {mustBeGreaterThan(NameValueArgs.stepSizePidParams,0)}
        NameValueArgs.stepSizePidTol (1,1) double {mustBePositive}
        NameValueArgs.stopOnBifurcation {validation.scalarLogical} % #ison:bifurcation
        NameValueArgs.stopOnCrossing {validation.scalarLogical} % #ison:bifurcation
        NameValueArgs.suspendContinuationErrorCounter (1,1) double {mustBeInteger} % #scalar#integer#positive#nonzero#larger:Opt.stepbackErrorCounter+1#neq:Opt.homotopyErrorCounter#neq:Opt.deflationErrorCounter
        NameValueArgs.targetTol (1,1) double {mustBePositive}
		NameValueArgs.timeLimit (1,1) double {mustBePositive,mustBeGreaterThan(NameValueArgs.timeLimit,0)}
        NameValueArgs.weightsAngleCustom (1,2) double {mustBePositive}
        NameValueArgs.weightsAngleChange (1,2) double {mustBePositive}
        NameValueArgs.weightsError (1,5) double {mustBePositive}
        NameValueArgs.weightsFayezioghani (1,2) double {mustBePositive}
        NameValueArgs.weightsMultiplicative (1,5) double {mustBePositive}
        NameValueArgs.weightsSzyszkowski (1,2) double {mustBePositive}
        NameValueArgs.weightsYoon (1,1) double {mustBePositive}
        NameValueArgs.xTarget (:,1) double
    end
    nameValueArgsCell = aux.struct2cellPreserveFieldnames(NameValueArgs);
    %
    %% initialize
    %
    warning on;
    [Opt,ds0,OptIsSet,func] = continuation.input(nameValueArgsCell,fun,var0,lStart,lEnd,ds0);
    [Opt,ds0,StepsizeOptions] = stepSize.initialize(Opt,var0,lStart,lEnd,ds0);
    if StepsizeOptions.rateOfContraction
        global solverStepsizes;
    end
    [Bifurcation,Counter,Do,Info,InfoOut,Initial,Is,Jacobian,Path,Plot,Plus,Remove,Solver,StepsizeInformation,Temp] = aux.initializeStructs(var0,lStart,lEnd,ds0,Opt,StepsizeOptions.rateOfContraction);
    clear('var0','lStart','lEnd','ds0');
    resCorr = continuation.corrector(fun,Opt);
    ds = Info.ds0;
    if ~Opt.jacobian
        aux.printLine(Opt,'Using finite differences to calculate Jacobian...\n');
    end
    aux.printLine(Opt,'Starting path continuation...\n');
    if Opt.dpaGammaVar
        paraName = 'g';
    else
        paraName = 'l';
    end
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
        % replace xPredictor with xTarget?:
        if OptIsSet.xTarget && norm(xPredictor-[Path.varAll(:,end);Path.lAll(end)])*1.1>norm(Opt.xTarget-[Path.varAll(:,end);Path.lAll(end)])
            xPredictor = Opt.xTarget;
            ds = norm(Opt.xTarget-[Path.varAll(:,end);Path.lAll(end)]);
            Do.convergeToxTarget = true;
            residual = @(x) aux.mergeResiduals(func,resCorr,x,[Path.varAll;Path.lAll],ds,Jacobian.last,Opt);
        end
        %
        %% solve
        %
        try
            dscale = aux.getDscale(Opt,Path);
            if sign(Path.lAll(end)-Opt.lTarget)*sign(lPredictor*(1+Opt.targetTol)-Opt.lTarget)<=0
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
        [Do,Info,Path,breakFunOut,Opt,Counter,ds] = aux.exitLoop(Do,Info,Initial,Is,Path,Opt,Counter,Bifurcation,ds,funSolution,Jacobian,breakFunOut);
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
    if Opt.dpa && OptIsSet.dpa && ~Opt.dpaGammaVar && Opt.bifurcation.parameterTrace
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
    aux.printLine(Opt,'--> time elapsed: %.3f s\n',toc(Info.t0));
    %
    %% output
    %
    InfoOut.numberOfSteps = Counter.step;
    InfoOut.numberOfInvalidPoints = Counter.loop - Counter.step;
    if OptIsSet.pathInfoFunction
        InfoOut.pathInfoValue = Path.pathInfoValue;
    end
    jacobianOut = Jacobian.last;
    exitflag = Info.exitflag;
    varAll = Path.varAll;
    lAll = Path.lAll;
    sAll = Path.sAll;
	%
end