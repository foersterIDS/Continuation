%% path continuation - bifurcation.check
%
%   Institute of Dynamics and Vibration Research
%   Leibniz University Hannover
%   26.05.2020 - Alwin Förster
%   02.07.2021 - Tido Kubatschek
%   25.05.2023 - Anna Lefken
%
function check(func,oih,resCorr)
    oih.bifurcation.flag = 0;
    lastJacobianRed = oih.path.getJacobianByName('last',1:oih.info.nv,1:oih.info.nv);
    lastJacobianRed = full(lastJacobianRed);
    bifFound=0;    % skip additional testfunction, when det(jac)=0
    if aux.ison(oih.opt.bifurcation) && oih.path.nAll>=3
        basicFoldDetection = (sign(diff(oih.path.lAll(end+[-1,0]))*diff(oih.path.lAll(end+[-2,-1])))<0);
        if ~isempty(oih.bifurcation.bif) && oih.bifurcation.bif(1,end)==oih.path.nAll-1
            basicFoldDetection = false;
        end
    else
        basicFoldDetection = false;
    end
    if oih.opt.bifurcation.mark
        %% mark bifurcations:
        signDetCurrentJacobianRed = sign(det(lastJacobianRed));
        if (signDetCurrentJacobianRed*oih.path.signDetJRedAll(end-1)<=0 && all(size(oih.path.getJacobianByName('previous'))==size(oih.path.getJacobianByName('last')))) || basicFoldDetection
            bifFound=1;
            bifType = (sign(det(oih.path.getJacobianByName('previous',1:oih.path.nVar,1:oih.path.nVar)))==sign(det(oih.path.getJacobianByName('last',1:oih.path.nVar,1:oih.path.nVar)))); % 1: fold bif.; 0: branch point bif; NaN: unknown
            oih.bifurcation.bif = [oih.bifurcation.bif,[oih.path.nAll;bifType]];
            oih.bifurcation.flag = 1;
            oih.bifurcation.scaling = [oih.bifurcation.scaling,1];
        end
    elseif oih.opt.bifurcation.determine || oih.opt.bifurcation.trace || oih.opt.bifurcation.parameterTrace
        %% determine bifurcation-points:
        jacStateTemp = oih.opt.jacobian;
        oih.opt.jacobian = false;
        [bifSolver,defaultBifSolverOutput] = continuation.solver(oih,0);
        oih.opt.jacobian = jacStateTemp;
        detSolverJacobianRed = det(lastJacobianRed);
        oih.bifurcation.scaling = [oih.bifurcation.scaling,1/detSolverJacobianRed];
        residualBif = @(x) bifurcation.residual(func,x,oih);
        signDetCurrentJacobianRed = sign(detSolverJacobianRed);
        if (signDetCurrentJacobianRed*oih.path.signDetJRedAll(end-1)<=0) || basicFoldDetection
            bifFound=1;
            %% find exact point:
            nds = 1000;
            indBif = oih.path.nAll;
            bifType = NaN;
            nv = numel(oih.path.varAll(:,1));
            
            if signDetCurrentJacobianRed*oih.path.signDetJRedAll(end-1)<=0
                dss = linspace(oih.path.sAll(end-1)-oih.path.sAll(end),0,nds);
                detLeft = det(oih.path.getJacobianByName('previous',1:nv,1:nv));
                detRight = det(oih.path.getJacobianByName('last',1:nv,1:nv));
                startDsp = detRight/(detRight-detLeft)*(dss(1)-dss(end));     % estimated bifurcation point
                indStart = find(dss>startDsp,1);                                % index of estimated bifurcation point
            else
                dss = linspace(oih.path.sAll(end-2)-oih.path.sAll(end),0,nds);
                dll = interp1(oih.path.sAll-oih.path.sAll(end),abs(oih.path.lAll-oih.path.lAll(end-1)),dss,'spline');
                [~,indStart] = min(dll);
            end
            
            indDss = zeros(nds,1);                                           % index vector for sorting dss
            indDss(1) = indStart;                                           % alternating left and right of indStart
            if indStart==1
                indDss = 1:1:nds;
            elseif indStart==nds
                indDss = nds:-1:1;
            elseif le(indStart,nds/2)
                indAlternatingLess = indStart-1:-1:1;
                indDss(2:2:length(indAlternatingLess)*2) = indAlternatingLess;
                indDss(3:2:length(indAlternatingLess)*2+1) = indStart+1:1:indStart+length(indAlternatingLess);
                indDss(length(indAlternatingLess)*2+2:end) = indStart+length(indAlternatingLess)+1:1:nds;
            else
                indAlternatingGreater = indStart+1:1:nds;
                indDss(2:2:length(indAlternatingGreater)*2) = indStart-1:-1:indStart-length(indAlternatingGreater);
                indDss(3:2:length(indAlternatingGreater)*2+1) = indAlternatingGreater;
                indDss(length(indAlternatingGreater)*2+2:nds) = indStart-length(indAlternatingGreater)-1:-1:1;
            end
            dss = dss(indDss);                                               % sort dss
            
            bifType = NaN;
            for ii=1:nds
                dsp = dss(ii);
                [varBifPredictor,lBifPredictor] = continuation.predictor(oih,dsp,lastJacobianRed,func,resCorr);
                dscale = aux.getDscale(oih);
                [xBif,funBif,bifSolverExitflag,bifSolverOutput,bifSolverJacobian] = bifSolver(residualBif,[varBifPredictor;lBifPredictor],dscale);
                if bifSolverExitflag>0
                    optAddPointArgs = {};
                    if oih.optIsSet.bifAdditionalTestfunction
                        optAddPointArgs{numel(optAddPointArgs)+1} = 'bifTestValue';
                        optAddPointArgs{numel(optAddPointArgs)+1} = oih.opt.bifAdditionalTestfunction(func,xBif,bifSolverJacobian,oih.path,oih.info);
                    end
                    if oih.optIsSet.pathInfoFunction
                        optAddPointArgs{numel(optAddPointArgs)+1} = 'pathInfoValue';
                        optAddPointArgs{numel(optAddPointArgs)+1} = oih.opt.pathInfoFunction(func,bifSolverJacobian,xBif(1:(end-1)),xBif(end));
                    end
                    if oih.stepsizeOptions.predictor
                        optAddPointArgs{numel(optAddPointArgs)+1} = 'predictor';
                        optAddPointArgs{numel(optAddPointArgs)+1} = [varBifPredictor;lBifPredictor];
                    end
                    if oih.stepsizeOptions.speedOfContinuation
                        optAddPointArgs{numel(optAddPointArgs)+1} = 'speedOfContinuation';
                        optAddPointArgs{numel(optAddPointArgs)+1} = 1;
                    end
                    oih.path.addPoint(xBif(1:end-1),xBif(end),bifSolverJacobian,oih,oih.path.nAll,optAddPointArgs{:});
                    % 
                    % get type of bifurcation
                    bifType = (sign(oih.path.detJv('name','previous'))==sign(oih.path.detJv('name','last'))); % 1: fold bif.; 0: branch point bif; NaN: unknown
                    %
                    % if bifType = 0 (branch point) calculate directions
                    % of paths by null() and save to oih.bifurcation.dirs cell array
                    if bifType == 0
                        nv = oih.path.nVar;
                        solverJacobianRed = bifSolverJacobian(1:nv,1:nv);
                        solverJacobianLam = bifSolverJacobian(1:nv,nv+1);
                        oih.bifurcation.jacobian = [oih.bifurcation.jacobian, {bifSolverJacobian}];
                        jacRedJacLam = [solverJacobianRed, solverJacobianLam];
                        %
                        nullMatrix = null(jacRedJacLam);
                        for kDirs = 1:width(nullMatrix)
                            bifDir = nullMatrix(:,kDirs);
                            if height(oih.bifurcation.dirs) == 1 && isempty(oih.bifurcation.dirs{1,2})
                                oih.bifurcation.dirs(1,:) = [{indBif}, {bifDir}];
                            else
                                oih.bifurcation.dirs(end+1,:) = [{indBif}, {bifDir}];
                            end
                        end
                    else
                        oih.bifurcation.jacobian = [oih.bifurcation.jacobian, {[]}];
                    end
                    break;
                end
            end
            oih.bifurcation.bif = [oih.bifurcation.bif,[indBif;bifType]];
            oih.bifurcation.flag = 1;
        end
    else
        %% error
        error('No such bifurcation-mode');
    end
    if oih.optIsSet.bifAdditionalTestfunction && bifFound==0
        if oih.opt.bifurcation.mark
            %% mark bifurcations:
            if oih.path.bifTestValue(end-1)*oih.path.bifTestValue(end)<=0
                bifFound=1;
                bifType = 2; % 2: additional; 1: fold bif.; 0: branch point bif; NaN: unknown
                oih.bifurcation.bif = [oih.bifurcation.bif,[oih.path.nAll;bifType]];
                oih.bifurcation.flag = 1;
                oih.bifurcation.scaling = [oih.bifurcation.scaling,1];
            end
        elseif oih.opt.bifurcation.determine || oih.opt.bifurcation.trace || oih.opt.bifurcation.parameterTrace
            %% determine bifurcation-points:
            jacStateTemp = oih.opt.jacobian;
            oih.opt.jacobian = false;
            [bifSolver,defaultBifSolverOutput] = continuation.solver(OptBif,0);
            oih.opt.jacobian = jacStateTemp;
            detSolverJacobianRed = det(lastJacobianRed);
            oih.bifurcation.scaling = [oih.bifurcation.scaling,1/detSolverJacobianRed];
            if oih.path.bifTestValue(end-1)*oih.path.bifTestValue(end)<=0
                bifFound=1;
                residualBif = @(x) bifurcation.residualAdditionalTestfunction(func,x,oih);
                %% find exact point:
                nds = 11;
                dss = linspace(oih.path.sAll(end-1)-oih.path.sAll(end),0,nds);
                dss = dss([6,7,5,8,4,9,3,10,2,11,1]);
                indBif = oih.path.nAll;
                bifType = NaN;
                for ii=1:nds
                    dsp = dss(ii);
                    [varBifPredictor,lBifPredictor] = continuation.predictor(oih,dsp,lastJacobianRed,func,resCorr);
                    dscale = aux.getDscale(oih);
                    [xBif,funBif,bifSolverExitflag,bifSolverOutput,bifSolverJacobian] = bifSolver(residualBif,[varBifPredictor;lBifPredictor],dscale);
                    if bifSolverExitflag>0
                        xBif=real(xBif);            % very small imaginary parts (probably caused by numeric differences)
                                                    % lead to errors in following code
                        oih.path.sAll = [oih.path.sAll(1:end-1),oih.path.sAll(end-1)+[norm(xBif-[oih.path.varAll(:,end-1);oih.path.lAll(end-1)]),norm(xBif-[oih.path.varAll(:,end-1);oih.path.lAll(end-1)])+norm([oih.path.varAll(:,end);oih.path.lAll(end)]-xBif)]];
                        oih.path.varAll = [oih.path.varAll(:,1:end-1),xBif(1:end-1),oih.path.varAll(:,end)];
                        oih.path.lAll = [oih.path.lAll(1:end-1),xBif(end),oih.path.lAll(end)];
                        % jacobian for xBif
                        Jacobianbif.solver=aux.numericJacobian(@(x)func(x(1:end-1),x(end)), xBif,'diffquot', oih.opt.diffquot); 
                        oih.path.bifTestValue=[oih.path.bifTestValue(1:end-1) oih.opt.bifAdditionalTestfunction(func,xBif,Jacobianbif,oih.path,oih.info),oih.path.bifTestValue(end)];
                        % 
                        % get type of bifurcation
                        bifType = 2; % 2: additional; 1: fold bif.; 0: branch point bif; NaN: unknown
                        break;
                    end
                end
                oih.bifurcation.bif = [oih.bifurcation.bif,[indBif;bifType]];
                oih.bifurcation.flag = 1;
            end
        else
            %% error
            error('No such bifurcation-mode');
        end
    end
end