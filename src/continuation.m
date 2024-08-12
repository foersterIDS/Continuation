%% path continuation - continuation
%
%   Institute of Dynamics and Vibration Research
%   Leibniz University Hannover
%   08.05.2020 - Alwin F�rster
%
%   [varAll,lAll,exitflag,Bifurcation,sAll,jacOut,breakFunOut,InfoOut] =...
%      continuation(fun,var0,lStart,lEnd,ds0,NameValueArgs)
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
% Continuation Copyright (C) 2024 Alwin F�rster
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
function [varAll,lAll,exitflag,bifStruct,sAll,jacobianOut,breakFunOut,infoOutStruct] = ...
    continuation(fun,var0,lStart,lEnd,ds0,NameValueArgs)
    %% arguments
    %
    arguments
        fun (1,1) function_handle
        var0 (:,1) double
        lStart (:,1) double
        lEnd (1,1) double
        ds0 (1,1) double {mustBeGreaterThan(ds0,0)}
        NameValueArgs.Opt (1,1) struct
        NameValueArgs.adaptCorrector (1,:) char {mustBeMember(NameValueArgs.adaptCorrector,{'basic','solve'})}
        NameValueArgs.alphaReverse (1,1) double {mustBeGreaterThan(NameValueArgs.alphaReverse,0)}
        NameValueArgs.alphaReverseAutoMode {validation.scalarLogical}
        NameValueArgs.approveManually {validation.scalarLogical} % #scalar#logical#ison:plot#ison:display|#scalar#double#ison:plot#ison:display
        NameValueArgs.bidirectional {validation.scalarLogical}
        NameValueArgs.bifurcation (1,:) char {mustBeMember(NameValueArgs.bifurcation,{'on','off','mark','determine','trace','parameterTrace'})}
        NameValueArgs.bifRandDir {validation.scalarLogical}
        NameValueArgs.bifAdditionalTestfunction (1,1) function_handle
        NameValueArgs.bifResidual (1,:) char {mustBeMember(NameValueArgs.bifResidual,{'determinant','luFactorization'})}
        NameValueArgs.breakFunction (1,1) function_handle
        NameValueArgs.checkJacobian {validation.scalarLogical}
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
        NameValueArgs.dsMin (1,1) double {mustBeGreaterThan(NameValueArgs.dsMin,0)} % #scalar#double#positive#smaller:norm(oih.opt.dsMax)
        NameValueArgs.dsTol (1,2) double {mustBePositive} % #array#positive#double#increasing#size:[1,2]
        NameValueArgs.dscale0 (:,1) double {mustBeGreaterThan(NameValueArgs.dscale0,0)} % #array#positive#nonzero#double#size:[numel(var0)+1,1]
        NameValueArgs.dscaleMin (:,1) double {mustBeGreaterThan(NameValueArgs.dscaleMin,0)} % #scalar#positive#nonzero#double|#array#positive#nonzero#double#size:[numel(var0)+1,1]
        NameValueArgs.enforceDsMax {validation.scalarLogical} % #ison:predictorSolver
        NameValueArgs.eventUserInput (1,:) cell
        NameValueArgs.g0 (1,1) double
        NameValueArgs.gTarget (1,1) double
        NameValueArgs.homotopy (1,:) char {mustBeMember(NameValueArgs.homotopy,{'on','off','f2','fix','fixnt','newton','squared'})}
        NameValueArgs.homotopyErrorCounter (1,1) double {mustBeInteger,mustBeGreaterThan(NameValueArgs.homotopyErrorCounter,0)} % #neq:oih.opt.deflationErrorCounter
        NameValueArgs.includeReverse {validation.scalarLogical}
        NameValueArgs.initialDeflationPoints (1,1) double
        NameValueArgs.jacobian {validation.scalarLogical}
        NameValueArgs.jacobianOut {mustBeMember(NameValueArgs.jacobianOut,{'basic','full'})}
        NameValueArgs.l0 (1,1) double
        NameValueArgs.lDirFunction (1,1) function_handle
        NameValueArgs.lMult0 (:,1) double
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
        NameValueArgs.pauseOnError {validation.scalarLogical}
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
        NameValueArgs.removeErrorCounter (1,1) double {mustBeInteger} % #scalar#integer#positive#nonzero#larger:oih.opt.stepbackErrorCounter+1#neq:oih.opt.homotopyErrorCounter#neq:oih.opt.deflationErrorCounter#neq:oih.opt.suspendContinuationErrorCounter|#scalar#integer#equals:0
        NameValueArgs.reverse  (1,:) char {mustBeMember(NameValueArgs.reverse,{'angle','jacobian'})}
        NameValueArgs.scaling (1,:) char {mustBeMember(NameValueArgs.scaling,{'dynamicdscale','staticdscale'})}
        NameValueArgs.solver (1,:) char {mustBeMember(NameValueArgs.solver,{'fsolve','lsqnonlin','newton'})}
        NameValueArgs.solverForce1it {validation.scalarLogical}
        NameValueArgs.solverMaxIterations (1,1) double {mustBeGreaterThan(NameValueArgs.solverMaxIterations,0),mustBeInteger}
        NameValueArgs.solverTol (1,1) double {mustBeGreaterThan(NameValueArgs.solverTol,0)}
        NameValueArgs.speedOfContinuation (1,1) double {mustBeGreaterThan(NameValueArgs.speedOfContinuation,0)}
        NameValueArgs.stepbackErrorCounter (1,1) double {mustBePositive,mustBeInteger,mustBeGreaterThan(NameValueArgs.stepbackErrorCounter,0)} % #neq:oih.opt.homotopyErrorCounter#neq:oih.opt.deflationErrorCounter
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
        NameValueArgs.suspendContinuationErrorCounter (1,1) double {mustBeInteger} % #scalar#integer#positive#nonzero#larger:oih.opt.stepbackErrorCounter+1#neq:oih.opt.homotopyErrorCounter#neq:oih.opt.deflationErrorCounter
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
    [oih,ds0,func] = continuation.input(nameValueArgsCell,fun,var0,lStart,lEnd,ds0);
    ds0 = stepSize.initialize(oih,var0,lStart,lEnd,ds0);
    if oih.stepsizeOptions.rateOfContraction
        global solverStepsizes;
    end
    oih.initializeStructsAndClasses(var0,lStart,lEnd,ds0);
    clear('var0','lStart','lEnd','ds0');
    resCorr = continuation.corrector(fun,oih);
    ds = oih.info.ds0;
    if ~oih.opt.jacobian
        aux.printLine(oih,'Using finite differences to calculate Jacobian...\n');
    end
    aux.printLine(oih,'Starting path continuation...\n');
    if oih.opt.dpaGammaVar
        paraName = 'g';
    else
        paraName = 'l';
    end
    speedOfContinuation = [];
    %
    %% find initial solution
    %
    residualInitial = @(v) aux.residualFixedValue(func,v,oih.opt.l0,oih);
    [varInitial,funInitial,initialExitflag,oih.solver.output,jacobianInitial] = oih.solver.main(residualInitial,oih.info.var0,oih.opt.dscale0(1:end-1));
    oih.solver.jacobian = jacobianInitial;
    aux.checkJacobian(residualInitial,funInitial,varInitial,oih);
    breakFunOut = [];
    event.eventObj = [];
    if initialExitflag>=0
        oih.do.continuation = true;
        oih.do.loop = true;
        dpaPoints = [];
        oih.solver.jacobian = [oih.solver.jacobian,aux.numericJacobian(@(x) func(x(1:oih.info.nv),x(oih.info.nv+1)),[varInitial;oih.info.lStart],'centralValue',funInitial,'derivativeDimensions',oih.info.nv+1,'diffquot',oih.opt.diffquot)];
        [~,breakFunOut] = oih.opt.breakFunction(funInitial,oih.solver.jacobian,varInitial,oih.info.lStart,breakFunOut);
        aux.printLine(oih,'Initial solution at %s = %.2e\n',paraName,oih.opt.l0);
        if oih.stepsizeOptions.rateOfContraction
            if size(solverStepsizes, 1) < 3
                oih.solver.output.rateOfContraction = oih.opt.optimalContractionRate;
            else
                oih.solver.output.rateOfContraction = solverStepsizes(3,2)/solverStepsizes(2,2);
            end
        end
        %% at initial point to path
        optArgsAddPoint = {};
        if oih.optIsSet.bifAdditionalTestfunction
            optArgsAddPoint = [optArgsAddPoint,'bifTestValue',1];
        end
        if oih.optIsSet.pathInfoFunction
            optArgsAddPoint = [optArgsAddPoint,'pathInfoValue',oih.opt.pathInfoFunction(func,oih.path,varInitial,oih.opt.l0)];
        end
        if oih.stepsizeOptions.predictor
            optArgsAddPoint = [optArgsAddPoint,'predictor',[oih.info.var0;oih.info.lStart]];
        end
        if oih.stepsizeOptions.speedOfContinuation
            optArgsAddPoint = [optArgsAddPoint,'speedOfContinuation',oih.opt.speedOfContinuation];
        end
        oih.path.addPointAtEnd(varInitial,oih.opt.l0,oih.solver.jacobian,oih,optArgsAddPoint{:});
        %% init. plot
        if aux.ison(oih.opt.plot)
            plot.livePlot(oih, oih.info.ds0, oih.info.ds0, oih.solver.output.iterations(end));
        end
    elseif oih.info.validJacobian
        oih.info.exitflag = -2;
        oih.do.continuation = false;
        aux.printLine(oih,'No initial solution found.\n');
    else
        aux.printLine(oih,'Provided Jacobian is corrupted.\n');
    end
    %
    %% continuation
    %
    while oih.do.continuation
        %% initialize loop
        %
        oih.counter.loop = oih.counter.loop+1;
        oih.is.currentJacobian = false;
        oih.counter.catchOld = oih.counter.catch;
        %
        %% residual
        %
        if oih.do.deflate
            try
                residual = @(x) aux.deflation(residual,xDeflation,x,oih.opt.jacobian);
            catch exceptionDeflation
                aux.printLine(oih,'---> deflation: catch!\n');
                if oih.opt.pauseOnError
                    aux.printLine(oih,['----> ',exceptionDeflation.message]);
                    input('Pause on error. Press [ENTER] to continue...');
                end
                oih.counter.catch = oih.counter.catch + 1;
            end
        elseif oih.do.suspend
            residual = @(v,lFix) aux.residualSuspendContinuation(func,v,lFix,oih);
        else
            residual = @(x) aux.mergeResiduals(func,resCorr,x,oih.path.xAll,ds,oih.path.getJacobianByName('last'),oih);
        end
        if oih.path.nAll==1 && ~isempty(oih.opt.initialDeflationPoints)
            for ii=1:numel(oih.opt.initialDeflationPoints(1,:))
                residual = @(x) aux.deflation(residual,oih.opt.initialDeflationPoints(:,ii),x,oih.opt.jacobian);
            end
        end
        %
        %% predictor
        %
        if oih.stepsizeOptions.speedOfContinuation
            timeNeeded = tic;
        end
        %
        try
            if oih.do.homotopy
                %% Homotopy
                xPredictor = homotopy.hContinuation(residual,[oih.path.varAll(:,end);oih.path.lAll(end)],oih);
                varPredictor = xPredictor(1:end-1);
                lPredictor = xPredictor(end);
            else
                %% calc. predictor
                [varPredictor,lPredictor,funPredictor,sPredictor,ds] = continuation.predictor(oih,ds,oih.path.getJacobianByName('last'),func,resCorr);
                xPredictor = [varPredictor;lPredictor];
            end
        catch exceptionPredictor
            [varPredictor,lPredictor,funPredictor,sPredictor,ds] = continuation.predictor(oih,ds,oih.path.getJacobianByName('last'),func,resCorr);
            xPredictor = [varPredictor;lPredictor];
            aux.printLine(oih,'---> predictor: catch!\n');
            if oih.opt.pauseOnError
                aux.printLine(oih,['----> ',exceptionPredictor.message]);
                input('Pause on error. Press [ENTER] to continue...');
            end
            oih.counter.catch = oih.counter.catch + 1;
        end
        % replace xPredictor with xTarget?:
        if oih.optIsSet.xTarget && norm(xPredictor-[oih.path.varAll(:,end);oih.path.lAll(end)])*1.1>norm(oih.opt.xTarget-[oih.path.varAll(:,end);oih.path.lAll(end)])
            xPredictor = oih.opt.xTarget;
            ds = norm(oih.opt.xTarget-[oih.path.varAll(:,end);oih.path.lAll(end)]);
            oih.do.convergeToxTarget = true;
            residual = @(x) aux.mergeResiduals(func,resCorr,x,oih.path.xAll,ds,oih.path.getJacobianByName('last'),oih);
        end
        %
        %% solve
        %
        try
            dscale = aux.getDscale(oih);
            if sign(oih.path.lAll(end)-oih.opt.lTarget)*sign(lPredictor*(1+oih.opt.targetTol)-oih.opt.lTarget)<=0
                %% try to converge to target
                residualTarget = @(v) aux.residualFixedValue(func,v,oih.opt.lTarget,oih);
                varPredictorCtt = (varPredictor - oih.path.varAll(:,end))*(abs(oih.opt.lTarget-oih.path.lAll(end))/abs(lPredictor-oih.path.lAll(end)))+oih.path.varAll(:,end);
                [varSolution,funSolution,oih.solver.exitflag,oih.solver.output,oih.solver.jacobian] = oih.solver.main(residualTarget,varPredictorCtt,dscale(1:end-1));
                xSolution = [varSolution;oih.opt.lTarget];
                oih.do.convergeToTarget = true;
                aux.checkJacobian(residualTarget,funSolution,varSolution,oih);
            else
                %% regular solver
                %            
                if oih.opt.solver.fsolve && oih.opt.solverForce1it
                    [xSolution,funSolution,oih.solver.exitflag,oih.solver.output,oih.solver.jacobian] = oih.solver.main(residual,xPredictor,dscale);
                    if oih.solver.output.iterations(end) < 1
                        % perturbate initial solution by tolerance of
                        % solver
                        pert = oih.opt.solverTol * ones(numel(xPredictor),1) / numel(xPredictor);
                        [xSolution,funSolution,oih.solver.exitflag,oih.solver.output,oih.solver.jacobian] = oih.solver.main(residual,xPredictor + pert,dscale);
                    end
                else
                    if oih.do.suspend
                        [vSolution,funSolution,oih.solver.exitflag,oih.solver.output,oih.solver.jacobian] = oih.solver.main(@(v) residual(v,xPredictor(end)),xPredictor(1:(end-1)),dscale(1:(end-1)));
                        xSolution = [vSolution;xPredictor(end)];
                    else
                        [xSolution,funSolution,oih.solver.exitflag,oih.solver.output,oih.solver.jacobian] = oih.solver.main(residual,xPredictor,dscale);
                    end
                end
                oih.do.convergeToTarget = false;
                aux.checkJacobian(residual,funSolution,xSolution,oih);
            end
            oih.is.currentJacobian = true;
        catch exceptionSolve
            xSolution = NaN(size(xPredictor));
            funSolution = inf(size(xPredictor));
            oih.solver.exitflag = -2;
            oih.solver.output = oih.solver.defaultOutput;
            oih.do.convergeToTarget = false;
            aux.printLine(oih,'---> solve: catch!\n');
            if oih.opt.pauseOnError
                aux.printLine(oih,['----> ',exceptionSolve.message]);
                input('Pause on error. Press [ENTER] to continue...');
            end
            oih.counter.catch = oih.counter.catch + 1;
        end
        %
        %% calc stepsize information
        %
        % measure speed
        %
        if oih.stepsizeOptions.speedOfContinuation
            timeNeeded = toc(timeNeeded);
            speedOfContinuation = ds/timeNeeded; 
        end
        %
        % measure rate of contraction
        if oih.stepsizeOptions.rateOfContraction
            if size(solverStepsizes, 1) < 3
                oih.solver.output.rateOfContraction = oih.opt.optimalContractionRate;
            else
                oih.solver.output.rateOfContraction = solverStepsizes(3,2)/solverStepsizes(2,2);
            end
        end
        %
        %% adaptive corrector
        %
        corrInfo = corrector.adapt(oih,func,xPredictor,dscale,oih.path.getJacobianByName('last'),ds);
        %
        %% check result
        %
        % check result:
        invPoiStr = aux.validateResult(xSolution,funSolution,oih,ds,funPredictor,sPredictor);
        % confirm result:
        xDeflation = aux.confirmResult(func,funSolution,xSolution,xPredictor,oih,speedOfContinuation);
        %
        %% Bifurcations
        %
        oih.bifurcation.flag = 0;
        if aux.ison(oih.opt.bifurcation) && oih.is.valid && ~oih.do.homotopy && oih.path.nAll>2
            if ~oih.is.currentJacobian
                %% get jacobian if not current
                oih.solver.jacobian = aux.getJacobian(func,oih.path.varAll(:,end),oih.path.lAll(end),oih);
            end
            bifurcation.check(func,oih,resCorr);
        end
        %
        %% DPA points
        %
        if oih.opt.dpa && oih.optIsSet.dpa && ~oih.opt.dpaGammaVar
            [dpaPoints] = dpa.checkResidual(fun,dpaPoints,oih);
        end
        %
        %% step size control
        %
        % save latest stepsize:
        dsim1 = ds;
        % adjust stepsize:
        [ds,event] = stepSize.control(ds,oih,event);
        %
        %% end loop
        %
        if oih.is.valid
            if ~isempty(oih.solver.output.iterations)
                aux.printLine(oih,'-----> continued at %s = %.4e\t|\tnew step size: ds = %.2e\t|\tloop counter = %d\t|\tstep = %d\t|\titerations = %d/%d\n',paraName,oih.path.lAll(end),ds,oih.counter.loop,oih.counter.step,oih.solver.output.iterations(end),oih.opt.nIterOpt);
            else
                aux.printLine(oih,'-----> continued at %s = %.4e\t|\tnew step size: ds = %.2e\t|\tloop counter = %d\t|\tstep = %d\t|\titerations = %d/%d\n',paraName,oih.path.lAll(end),ds,oih.counter.loop,oih.counter.step,[],oih.opt.nIterOpt);
            end
        else
            if ~isempty(oih.solver.output.iterations)
                aux.printLine(oih,'-----> invalid point %s |\tnew step size: ds = %.2e\t|\tloop counter = %d\t|\tstep = %d\t|\titerations = %d/%d\n',invPoiStr,ds,oih.counter.loop,oih.counter.step,oih.solver.output.iterations(end),oih.opt.nIterOpt);
            else
                aux.printLine(oih,'-----> invalid point %s |\tnew step size: ds = %.2e\t|\tloop counter = %d\t|\tstep = %d\t|\titerations = %d/%d\n',invPoiStr,ds,oih.counter.loop,oih.counter.step,[],oih.opt.nIterOpt);
            end
        end
        [breakFunOut,ds] = aux.exitLoop(oih,ds,funSolution,breakFunOut);
        if oih.do.changeCorrector
            oih.opt = aux.seton(oih.opt,'corrector',corrInfo);
            resCorr = continuation.corrector(fun,oih);
            oih.do.changeCorrector = false;
        end
        oih.updateOpt();
        %
        %% live plot
        %
        if aux.ison(oih.opt.plot) && oih.is.valid
            try
                plot.livePlot(oih, ds, dsim1, oih.solver.output.iterations(end), funPredictor, sPredictor, dpaPoints);
            catch exceptionPlot
                aux.printLine(oih,'--> The plot update has failed.\n');
                if oih.opt.pauseOnError
                    aux.printLine(oih,['----> ',exceptionPlot.message]);
                    input('Pause on error. Press [ENTER] to continue...');
                end
            end
        end
        %
        %% catch counter
        %
        if oih.counter.catchOld == oih.counter.catch
            oih.counter.catch = 0;
        end
        %
    end
    %
    %% bifurcation tracing
    %
    if oih.opt.bifurcation.trace
        try
            delete(oih.plot.plCurr);
            bifurcation.trace(oih,func,resCorr);
        catch exceptionBifurcationTrace
            aux.printLine(oih,'--> Failed to trace bifurcations.\n');
            if oih.opt.pauseOnError
                aux.printLine(oih,['---> ',exceptionBifurcationTrace.message]);
                input('Pause on error. Press [ENTER] to continue...');
            end
        end
    elseif oih.opt.bifurcation.parameterTrace
        try
            delete(oih.plot.plCurr);
            bifurcation.parameterTrace(oih,fun);
        catch exceptionBifurcationParameterTrace
            aux.printLine(oih,'--> Failed to trace bifurcation parameter.\n');
            if oih.opt.pauseOnError
                aux.printLine(oih,['---> ',exceptionBifurcationParameterTrace.message]);
                input('Pause on error. Press [ENTER] to continue...');
            end
        end
    end
    %
    %% DPA
    %
    if oih.opt.dpa && oih.optIsSet.dpa && ~oih.opt.dpaGammaVar && oih.opt.bifurcation.parameterTrace
        Path = dpa.trace(fun,dpaPoints,oih);
    end
    %
    %% live plot finalization
    %
    if aux.ison(oih.opt.plot) && initialExitflag>0
        try
            oih.info.finalSolutionPoint = true;
            plot.livePlot(oih, ds, dsim1, oih.solver.output.iterations(end), funPredictor, sPredictor, dpaPoints);
            if isfield(oih.plot,'plCurr')
                delete(oih.plot.plCurr);
            end
        catch exceptionPlotFinal
            aux.printLine(oih,'--> The plot update has failed.\n');
            if oih.opt.pauseOnError
                aux.printLine(oih,['---> ',exceptionPlotFinal.message]);
                input('Pause on error. Press [ENTER] to continue...');
            end
        end
        drawnow;
    end
    %
    %% final disp
    %
    aux.printLine(oih,[oih.info.exitMsg,'\n']);
    aux.printLine(oih,'--> time elapsed: %.3f s\n',toc(oih.info.t0));
    %
    %% output
    %
    oih.infoOut.numberOfSteps = oih.counter.step;
    oih.infoOut.numberOfInvalidPoints = oih.counter.loop - oih.counter.step;
    if oih.optIsSet.pathInfoFunction
        oih.infoOut.pathInfoValue = oih.path.pathInfoValue;
    end
    if oih.opt.jacobianOut.basic
        jacobianOut = oih.path.getJacobianByName('last');
    elseif oih.opt.jacobianOut.full
        jacobianOut = oih.path.JAll;
    end
    exitflag = oih.info.exitflag;
    varAll = oih.path.varAll;
    lAll = oih.path.lAll;
    sAll = oih.path.sAll;
    bifStruct = oih.bifurcation;
    infoOutStruct = oih.infoOut;
	%
end