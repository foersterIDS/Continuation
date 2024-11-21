%% path continuation - aux.confirmResult
%
%   Institute of Dynamics and Vibration Research
%   Leibniz University Hannover
%   28.03.2022 - Alwin Förster
%   12.08.2022 - Anna Lefken
%
function [xDeflation] = confirmResult(func,funSolution,xSolution,xPredictor,oih,speedOfContinuation)
    xDeflation = xSolution;
    if oih.is.valid
        %% valid result
        % Add point
        optAddPointArgs = {};
        if oih.optIsSet.bifAdditionalTestfunction
            optAddPointArgs{numel(optAddPointArgs)+1} = 'bifTestValue';
            optAddPointArgs{numel(optAddPointArgs)+1} = oih.opt.bifAdditionalTestfunction(func,xSolution,oih.path,oih.info);
        end
        if oih.optIsSet.pathInfoFunction
            optAddPointArgs{numel(optAddPointArgs)+1} = 'pathInfoValue';
            optAddPointArgs{numel(optAddPointArgs)+1} = oih.opt.pathInfoFunction(func,oih.solver.jacobian,xSolution(1:(end-1)),xSolution(end));
        end
        if oih.stepsizeOptions.predictor
            optAddPointArgs{numel(optAddPointArgs)+1} = 'predictor';
            optAddPointArgs{numel(optAddPointArgs)+1} = xPredictor;
        end
        if oih.stepsizeOptions.speedOfContinuation
            optAddPointArgs{numel(optAddPointArgs)+1} = 'speedOfContinuation';
            optAddPointArgs{numel(optAddPointArgs)+1} = speedOfContinuation;
        end
        oih.path.addPointAtEnd(xSolution(1:(end-1)),xSolution(end),oih.solver.jacobian,oih,optAddPointArgs{:});
        % set structs
        if ~oih.path.stepBackStatus
            oih.counter.validStepback = 0;
        else
            oih.counter.validStepback = oih.counter.validStepback+1;
            oih.path.toggleStepback();
        end
        oih.do.deflate = false;
        oih.do.homotopy = false;
        oih.do.stepback = false;
        oih.do.suspend = false;
        oih.counter.error = 0;
        oih.counter.step = oih.counter.step + 1;
        if oih.do.remove
            if oih.path.sAll(end)>(oih.remove.s+oih.remove.ds)
                oih.opt.dsMax = oih.initial.dsMax;
                oih.do.remove = false;
                oih.counter.remove = 0;
            elseif oih.counter.step>5+oih.remove.step
                oih.opt.dsMax = oih.initial.dsMax;
                oih.do.remove = false;
            end
        end
    else
        %% invalid result
        oih.counter.error = oih.counter.error+1;
        if oih.opt.deflation && ~oih.do.deflate && ~isnan(sum(xSolution(:,end))) && oih.is.reverse && oih.counter.error>=oih.opt.deflationErrorCounter
            oih.do.deflate = true;
            xDeflation = xSolution;
        else
            oih.do.deflate = false;
        end
        if ((oih.counter.error==oih.opt.stepbackErrorCounter) || oih.do.stepbackManually) && (oih.path.nAll>1)
            %% stepback
            oih.path.toggleStepback();
            if oih.counter.validStepback<oih.opt.stepbackErrorCounter && ~oih.optIsSet.lMult0
                oih.do.stepback = true;
            else
                oih.do.stepback = false;
            end
        elseif (oih.counter.error==oih.opt.stepbackErrorCounter+1) && (oih.path.nAll>1)
            %% undo stepback
            oih.path.toggleStepback();
            oih.do.stepback = false;
            oih.do.suspend = false;
        elseif (oih.counter.error==oih.opt.suspendContinuationErrorCounter) && (oih.path.nAll>1)
            %% suspend
            oih.path.suspend();
            oih.do.stepback = false;
            oih.do.suspend = true;
        elseif logical(oih.opt.removeErrorCounter) && ((oih.counter.error==oih.opt.removeErrorCounter) && (oih.path.nAll>1))
            %% remove
            nPath = oih.path.nAll;
            oih.remove.s = oih.path.sAll(nPath);
            nRmv = min([2*oih.opt.removeErrorCounter,nPath-1]);
            oih.opt.dsMax = max([mean(diff(oih.path.sAll(nPath+((-ceil(nRmv/2)+1):0))))*0.75,oih.opt.dsMin,oih.opt.dsMin*10]);
            oih.path.remove(nPath+((-nRmv+1):0));
            oih.remove.ds = oih.remove.s-oih.path.sAll(end);
            dscaleRmv = aux.getDscale(oih);
            lRmv = oih.path.lAll(end);
            residualFixedValueRmv = @(v) aux.residualFixedValue(func,v,lRmv,oih);
            [varRmv,funRmv,~,~,rmvJacobian] = oih.solver.main(residualFixedValueRmv,oih.path.varAll(:,end),dscaleRmv(1:end-1));
            if ~isempty(oih.bifurcation.bif) && numel(oih.bifurcation.bif(1,:))>0
                nBifsRmv = sum(sum(oih.bifurcation.bif(1,:)'==(nPath+((-nRmv+1):0))));
                if nBifsRmv>0
                    oih.bifurcation.bif(:,end+((1-nBifsRmv):0)) = [];
                end
            end
            oih.path.remove(oih.path.nAll);
            xRmv = [varRmv;lRmv];
            optAddPointArgs = {};
            if oih.optIsSet.bifAdditionalTestfunction
                optAddPointArgs{numel(optAddPointArgs)+1} = 'bifTestValue';
                optAddPointArgs{numel(optAddPointArgs)+1} = oih.opt.bifAdditionalTestfunction(func,xRmv,rmvJacobian,oih.path,oih.info);
            end
            if oih.optIsSet.pathInfoFunction
                optAddPointArgs{numel(optAddPointArgs)+1} = 'pathInfoValue';
                optAddPointArgs{numel(optAddPointArgs)+1} = oih.opt.pathInfoFunction(func,rmvJacobian,xRmv(1:(end-1)),xRmv(end));
            end
            if oih.stepsizeOptions.predictor
                optAddPointArgs{numel(optAddPointArgs)+1} = 'predictor';
                optAddPointArgs{numel(optAddPointArgs)+1} = xRmv;
            end
            if oih.stepsizeOptions.speedOfContinuation
                optAddPointArgs{numel(optAddPointArgs)+1} = 'speedOfContinuation';
                optAddPointArgs{numel(optAddPointArgs)+1} = speedOfContinuation;
            end
            rmvJacobian = [rmvJacobian,aux.numericJacobian(@(x) func(x(1:oih.info.nv),x(oih.info.nv+1)),xRmv,'centralValue',funRmv,'derivativeDimensions',oih.info.nv+1,'diffquot',oih.opt.diffquot,'diffStep',oih.opt.diffStep)];
            oih.path.addPointAtEnd(varRmv,lRmv,rmvJacobian,oih,optAddPointArgs{:});
            oih.do.stepback = false;
            oih.do.suspend = false;
            oih.do.remove = true;
            oih.remove.step = oih.counter.step;
            oih.counter.remove = oih.counter.remove+1;
        else
            %% else
            oih.do.stepback = false;
            oih.do.suspend = false;
            oih.do.remove = false;
        end
        if oih.opt.includeReverse && oih.is.reverse && oih.solver.exitflag>0
            aux.includeReverse(xSolution,oih);
        end
        if aux.ison(oih.opt.homotopy) && oih.counter.error>=oih.opt.homotopyErrorCounter
            oih.do.homotopy = true;
        else
            oih.do.homotopy = false;
        end
        if oih.is.catch
            oih.counter.catch = oih.counter.catch + 1;
        end
    end
end