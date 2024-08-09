%% path continuation - aux.confirmResult
%
%   Institute of Dynamics and Vibration Research
%   Leibniz University Hannover
%   28.03.2022 - Alwin Förster
%   12.08.2022 - Anna Lefken
%
function [xDeflation,Bifurcation,Counter,Do,Info,Initial,Is,Remove,Solver,StepsizeOptions,Opt] = ...
    confirmResult(func,funSolution,xSolution,xPredictor,Bifurcation,Counter,Do,Info,Initial,Is,Path,Remove,Solver,speedOfContinuation,StepsizeOptions,Opt,OptIsSet)
    xDeflation = xSolution;
    if Is.valid
        %% valid result
        % Add point
        optAddPointArgs = {};
        if OptIsSet.bifAdditionalTestfunction
            optAddPointArgs{numel(optAddPointArgs)+1} = 'bifTestValue';
            optAddPointArgs{numel(optAddPointArgs)+1} = Opt.bifAdditionalTestfunction(func,xSolution,Path,Info);
        end
        if OptIsSet.pathInfoFunction
            optAddPointArgs{numel(optAddPointArgs)+1} = 'pathInfoValue';
            optAddPointArgs{numel(optAddPointArgs)+1} = Opt.pathInfoFunction(func,Solver.jacobian,xSolution(1:(end-1)),xSolution(end));
        end
        if StepsizeOptions.predictor
            optAddPointArgs{numel(optAddPointArgs)+1} = 'predictor';
            optAddPointArgs{numel(optAddPointArgs)+1} = xPredictor;
        end
        if StepsizeOptions.speedOfContinuation
            optAddPointArgs{numel(optAddPointArgs)+1} = 'speedOfContinuation';
            optAddPointArgs{numel(optAddPointArgs)+1} = speedOfContinuation;
        end
        Path.addPointAtEnd(xSolution(1:(end-1)),xSolution(end),Solver.jacobian,Solver,optAddPointArgs{:});
        % set structs
        if ~Path.stepBackStatus
            Counter.validStepback = 0;
        else
            Counter.validStepback = Counter.validStepback+1;
            Path.toggleStepback();
        end
        Do.deflate = false;
        Do.homotopy = false;
        Do.stepback = false;
        Do.suspend = false;
        Counter.error = 0;
        Counter.step = Counter.step + 1;
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
        if ((Counter.error==Opt.stepbackErrorCounter) || Do.stepbackManually) && (Path.nAll>1)
            %% stepback
            Path.toggleStepback();
            if Counter.validStepback<Opt.stepbackErrorCounter
                Do.stepback = true;
            else
                Do.stepback = false;
            end
        elseif (Counter.error==Opt.stepbackErrorCounter+1) && (Path.nAll>1)
            %% undo stepback
            Path.toggleStepback();
            Do.stepback = false;
            Do.suspend = false;
        elseif (Counter.error==Opt.suspendContinuationErrorCounter) && (Path.nAll>1)
            %% suspend
            Path.suspend();
            Do.stepback = false;
            Do.suspend = true;
        elseif logical(Opt.removeErrorCounter) && ((Counter.error==Opt.removeErrorCounter) && (Path.nAll>1))
            %% remove
            nPath = Path.nAll;
            Remove.s = Path.sAll(nPath);
            nRmv = min([2*Opt.removeErrorCounter,nPath-1]);
            Opt.dsMax = max([mean(diff(Path.sAll(nPath+((-ceil(nRmv/2)+1):0))))*0.75,Opt.dsMin,Opt.dsMin*10]);
            Path.remove(nPath+((-nRmv+1):0));
            Remove.ds = Remove.s-Path.sAll(end);
            dscaleRmv = aux.getDscale(Opt,Path);
            lRmv = Path.lAll(end);
            residualFixedValueRmv = @(v) aux.residualFixedValue(func,v,lRmv,Opt);
            [varRmv,funRmv,~,~,rmvJacobian] = Solver.main(residualFixedValueRmv,Path.varAll(:,end),dscaleRmv(1:end-1));
            if ~isempty(Bifurcation.bif) && numel(Bifurcation.bif(1,:))>0
                nBifsRmv = sum(sum(Bifurcation.bif(1,:)'==(nPath+((-nRmv+1):0))));
                if nBifsRmv>0
                    Bifurcation.bif(:,end+((1-nBifsRmv):0)) = [];
                end
            end
            Path.remove(Path.nAll);
            xRmv = [varRmv;lRmv];
            optAddPointArgs = {};
            if OptIsSet.bifAdditionalTestfunction
                optAddPointArgs{numel(optAddPointArgs)+1} = 'bifTestValue';
                optAddPointArgs{numel(optAddPointArgs)+1} = Opt.bifAdditionalTestfunction(func,xRmv,rmvJacobian,Path,Info);
            end
            if OptIsSet.pathInfoFunction
                optAddPointArgs{numel(optAddPointArgs)+1} = 'pathInfoValue';
                optAddPointArgs{numel(optAddPointArgs)+1} = Opt.pathInfoFunction(func,rmvJacobian,xRmv(1:(end-1)),xRmv(end));
            end
            rmvJacobian = [rmvJacobian,aux.numericJacobian(@(x) func(x(1:Info.nv),x(Info.nv+1)),[Path.varAll(:,end);Path.lAll(end)],'centralValue',funRmv,'derivativeDimensions',Info.nv+1,'diffquot',Opt.diffquot)];
            Path.addPointAtEnd(varRmv,lRmv,rmvJacobian,optAddPointArgs{:});
            Do.stepback = false;
            Do.suspend = false;
            Do.remove = true;
            Counter.remove = Counter.remove+1;
        else
            %% else
            Do.stepback = false;
            Do.suspend = false;
            Do.remove = false;
        end
        if Opt.includeReverse && Is.reverse && Solver.exitflag>0
            aux.includeReverse(xSolution,Path,Solver);
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