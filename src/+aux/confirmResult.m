%% path continuation - aux.confirmResult
%
%   Institute of Dynamics and Vibration Research
%   Leibniz University Hannover
%   28.03.2022 - Alwin Förster
%   12.08.2022 - Anna Lefken
%
function [xDeflation,Bifurcation,Counter,Do,Info,Initial,Is,Jacobian,Path,Plus,Remove,Solver,StepsizeInformation,StepsizeOptions,Temp,Opt] = ...
    confirmResult(func,xSolution,xPredictor,Bifurcation,Counter,Do,Info,Initial,Is,Jacobian,Path,Plus,Remove,Solver,StepsizeInformation,StepsizeOptions,Temp,Opt,OptIsSet)
    xDeflation = xSolution;
    if Is.valid
        %% valid result
        if isempty(Plus.x)
            Path.varAll = [Path.varAll,xSolution(1:end-1)];
            Path.lAll = [Path.lAll,xSolution(end)];
            Path.sAll = [Path.sAll,Path.sAll(end)+norm(xSolution-[Path.varAll(:,end-1);Path.lAll(end-1)])];
            if OptIsSet.bifAdditionalTestfunction
                Path.biftestValue=[Path.biftestValue Opt.bifAdditionalTestfunction(func,xSolution,Jacobian,Path,Info)];
            end
            % update stepsize information
            % measure speed
            %
            if StepsizeOptions.speedOfContinuation
                if ~isempty(Temp.speedOfContinuation)
                    Temp.speedOfContinuation = [Temp.speedOfContinuation, StepsizeInformation.speedOfContinuation];
                else
                    Temp.speedOfContinuation = StepsizeInformation.speedOfContinuation;
                end
            end
            %
            % measure rate of contraction
            if StepsizeOptions.rateOfContraction
                if ~isempty(Temp.rateOfContraction)
                    Temp.rateOfContraction = [Temp.rateOfContraction, StepsizeInformation.rateOfContraction];
                else
                    Temp.rateOfContraction = StepsizeInformation.rateOfContraction;
                end
            end
            %
            if StepsizeOptions.predictor
                if ~isempty(Temp.predictor)
                    Temp.predictor = [Temp.predictor, xPredictor];
                end
            end
            Jacobian.previous = Jacobian.last;
            Jacobian.last = Jacobian.solver;
            Counter.validStepback = 0;
        else
            Path.varAll = [Path.varAll,xSolution(1:end-1),Plus.x(1:end-1)];
            Path.lAll = [Path.lAll,xSolution(end),Plus.x(end)];
            Path.sAll = [Path.sAll,Path.sAll(end)+norm(xSolution-[Path.varAll(:,end-2);Path.lAll(end-2)])*[1,1]+norm(Plus.x-xSolution)*[0,1]];
            if OptIsSet.bifAdditionalTestfunction
                Path.biftestValue=[Path.biftestValue,Opt.bifAdditionalTestfunction(func,xSolution,Jacobian,Path,Info),Plus.biftestValue];
            end
            
            if StepsizeOptions.predictor
                Temp.predictor = [Temp.predictor, xPredictor, Plus.xPredictor];
                Plus.xPredictor = [];
            end
            if StepsizeOptions.speedOfContinuation
                Temp.speedOfContinuation = [Temp.speedOfContinuation, StepsizeInformation.speedOfContinuation, Plus.speedOfContinuation];
                Plus.speedOfContinuation = [];
            end
            if StepsizeOptions.rateOfContraction
                Temp.rateOfContraction = [Temp.rateOfContraction, StepsizeInformation.rateOfContraction, Plus.rateOfContraction];
                Plus.rateOfContraction = [];
            end
            Plus.x = [];
            Jacobian.previous = Jacobian.solver;
            Jacobian.last = Plus.jacobian;
            Plus.jacobian = [];
            Counter.validStepback = Counter.validStepback+1;
        end
        Do.deflate = false;
        Do.homotopy = false;
        Do.stepback = false;
        Do.suspend = false;
        Counter.error = 0;
        Counter.step = Counter.step + 1;
        if Do.convergeToTarget
            Jacobian.last = aux.getJacobian(func,Path.varAll(:,end),Path.lAll(end),Opt);
        end
        if Do.remove
            if Path.sAll(end)>(Remove.s+Remove.ds)
                Opt.dsMax = Initial.dsMax;
                Do.remove = false;
                Counter.remove = 0;
            end
        end
    else
        %% invalid result
        Counter.error = Counter.error+1;
        if Opt.deflation && ~Do.deflate && ~isnan(sum(xSolution(:,end))) && Is.reverse && Counter.error>=Opt.deflationErrorCounter
            Do.deflate = true;
            xDeflation = xSolution;
        else
            Do.deflate = false;
        end
        if ((Counter.error==Opt.stepbackErrorCounter) || Do.stepbackManually) && (numel(Path.lAll)>1)
            %% stepback
            if Counter.validStepback<Opt.stepbackErrorCounter
                Do.stepback = true;
                Plus.x = [Path.varAll(:,end);Path.lAll(end)];
                if OptIsSet.bifAdditionalTestfunction
                    Plus.biftestValue=Opt.bifAdditionalTestfunction(func,xSolution,Jacobian,Path,Info);
                end
                if StepsizeOptions.predictor
                    Plus.xPredictor = xPredictor;
                end
                if StepsizeOptions.speedOfContinuation
                    Plus.speedOfContinuation = StepsizeInformation.speedOfContinuation;
                end
                if StepsizeOptions.rateOfContraction
                    Plus.rateOfContraction = StepsizeInformation.rateOfContraction;
                end
            else
                Do.stepback = false;
            end
            Path.varAll(:,end) = [];
            Path.lAll(end) = [];
            Path.sAll(end) = [];
            if OptIsSet.bifAdditionalTestfunction
                Path.biftestValue(end) = [];
            end
            if StepsizeOptions.predictor
                if ~isempty(Temp.predictor)
                    Temp.predictor(:,end) = [];
                else
                    Temp.predictor = [];
                end
            end
            if StepsizeOptions.speedOfContinuation && ~isempty(Temp.speedOfContinuation)
                Temp.speedOfContinuation(end) = [];
            end
            if StepsizeOptions.rateOfContraction && ~isempty(Temp.rateOfContraction)
                Temp.rateOfContraction(end) = [];
            end
            Plus.jacobian = Jacobian.last;
            Jacobian.last = Jacobian.previous;
        elseif (Counter.error==Opt.stepbackErrorCounter+1) && (numel(Path.lAll)>1)
            %% undo stepback
            if ~isempty(Plus.x)
                Path.varAll = [Path.varAll,Plus.x(1:end-1)];
                Path.lAll = [Path.lAll,Plus.x(end)];
                Path.sAll = [Path.sAll,Path.sAll(end)+norm([Path.varAll(:,end);Path.lAll(end)]-[Path.varAll(:,end-1);Path.lAll(end-1)])];
                if OptIsSet.bifAdditionalTestfunction
                    Path.biftestValue = [Path.biftestValue,Plus.biftestValue];
                end
                if StepsizeOptions.predictor
                    Temp.predictor = [Temp.predictor, Plus.xPredictor];
                end
                if StepsizeOptions.speedOfContinuation
                    Temp.speedOfContinuation = [Temp.speedOfContinuation, Plus.speedOfContinuation];
                end
                if StepsizeOptions.rateOfContraction
                    Temp.rateOfContraction = [Temp.rateOfContraction, Plus.rateOfContraction];
                end
                Plus = aux.clearStruct(Plus);
            end
            Do.stepback = false;
            Do.suspend = false;
        elseif (Counter.error==Opt.suspendContinuationErrorCounter) && (numel(Path.lAll)>1)
            %% suspend
            Plus = aux.clearStruct(Plus);
            Do.stepback = false;
            Do.suspend = true;
        elseif logical(Opt.removeErrorCounter) && ((Counter.error==Opt.removeErrorCounter) && (numel(Path.lAll)>1))
            %% remove
            nPath = numel(Path.lAll);
            Remove.s = Path.sAll(nPath);
            nRmv = min([2*Opt.removeErrorCounter,nPath-1]);
            Opt.dsMax = max([mean(diff(Path.sAll(nPath+((-ceil(nRmv/2)+1):0))))*0.75,Opt.dsMin,Opt.dsMin*10]);
            Path.varAll(:,nPath+((-nRmv+1):0)) = [];
            Path.lAll(nPath+((-nRmv+1):0)) = [];
            Path.sAll(nPath+((-nRmv+1):0)) = [];
            if OptIsSet.bifAdditionalTestfunction
                Path.biftestValue(nPath+((-nRmv+1):0)) = [];
            end
            Remove.ds = Remove.s-Path.sAll(end);
            dscaleRmv = aux.getDscale(Opt,Path);
            residualFixedValueRmv = @(v) aux.residualFixedValue(func,v,Path.lAll(end),Opt);
            [varRmv,funRmv,~,~,rmvJacobian] = Solver.main(residualFixedValueRmv,Path.varAll(:,end),dscaleRmv(1:end-1));
            if ~isempty(Bifurcation.bif) && numel(Bifurcation.bif(1,:))>0
                nBifsRmv = sum(sum(Bifurcation.bif(1,:)'==(nPath+((-nRmv+1):0))));
                if nBifsRmv>0
                    Bifurcation.bif(:,end+((1-nBifsRmv):0)) = [];
                    Jacobian.signDetRed = Jacobian.signDetRed*(-1)^(nBifsRmv);
                end
            end
            Path.varAll(:,end) = varRmv;
            Jacobian.last = [rmvJacobian,aux.numericJacobian(@(x) func(x(1:Info.nv),x(Info.nv+1)),[Path.varAll(:,end);Path.lAll(end)],'centralValue',funRmv,'derivativeDimensions',Info.nv+1,'diffquot',Opt.diffquot)];
            if StepsizeOptions.predictor
                if nRmv >= 3
                    Temp.predictor = xPredictor;
                else
                    Temp.predictor = Temp.predictor(:,1:end-nRmv);
                end
            end
            if StepsizeOptions.speedOfContinuation
                if nRmv >= 3
                    Temp.speedOfContinuation = [];
                else
                    Temp.speedOfContinuation = Temp.speedOfContinuation(1:end-nRmv);
                end
            end
            if StepsizeOptions.rateOfContraction
                if nRmv >= 3
                    Temp.rateOfContraction = [];
                else
                    Temp.rateOfContraction = Temp.rateOfContraction(1:end-nRmv);
                end
            end
            Plus = aux.clearStruct(Plus);
            Do.stepback = false;
            Do.suspend = false;
            Do.remove = true;
            Counter.remove = Counter.remove+1;
        else
            %% else
            Plus = aux.clearStruct(Plus);
            Do.stepback = false;
            Do.suspend = false;
            Do.remove = false;
        end
        if Opt.includeReverse && Is.reverse && Solver.exitflag>0
            Path = aux.includeReverse(xSolution,Path);
        end
        if aux.ison(Opt.homotopy) && Counter.error>=Opt.homotopyErrorCounter
            Do.homotopy = true;
        else
            Do.homotopy = false;
        end
        if Is.catch
            Counter.catch = Counter.catch + 1;
        end
    end
end