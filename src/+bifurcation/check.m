%% path continuation - bifurcation.check
%
%   Institute of Dynamics and Vibration Research
%   Leibniz University Hannover
%   26.05.2020 - Alwin Förster
%   02.07.2021 - Tido Kubatschek
%   25.05.2023 - Anna Lefken
%
function [Bifurcation] = check(func,Path,Bifurcation,Info,resCorr,Solver,Opt,OptIsSet)
    Bifurcation.flag = 0;
    lastJacobianRed = Path.getJacobianByName('last',1:Info.nv,1:Info.nv);
    lastJacobianRed = full(lastJacobianRed);
    bifFound=0;    % skip additional testfunction, when det(jac)=0
    if aux.ison(Opt.bifurcation) && Path.nAll>=3
        basicFoldDetection = (sign(diff(Path.lAll(end+[-1,0]))*diff(Path.lAll(end+[-2,-1])))<0);
        if ~isempty(Bifurcation.bif) && Bifurcation.bif(1,end)==Path.nAll-1
            basicFoldDetection = false;
        end
    else
        basicFoldDetection = false;
    end
    if Opt.bifurcation.mark
        %% mark bifurcations:
        signDetCurrentJacobianRed = sign(det(lastJacobianRed));
        if (signDetCurrentJacobianRed*Path.signDetJRedAll(end-1)<=0 && all(size(Path.getJacobianByName('previous'))==size(Path.getJacobianByName('last')))) || basicFoldDetection
            bifFound=1;
            bifType = (sign(det(Path.getJacobianByName('previous',1:Path.nVar,1:Path.nVar)))==sign(det(Path.getJacobianByName('last',1:Path.nVar,1:Path.nVar)))); % 1: fold bif.; 0: branch point bif; NaN: unknown
            Bifurcation.bif = [Bifurcation.bif,[Path.nAll;bifType]];
            Bifurcation.flag = 1;
            Bifurcation.scaling = [Bifurcation.scaling,1];
        end
    elseif Opt.bifurcation.determine || Opt.bifurcation.trace || Opt.bifurcation.parameterTrace
        %% determine bifurcation-points:
        OptBif = Opt;
        OptBif.jacobian = false;
        [bifSolver,defaultBifSolverOutput] = continuation.solver(OptBif,0);
        detSolverJacobianRed = det(lastJacobianRed);
        Bifurcation.scaling = [Bifurcation.scaling,1/detSolverJacobianRed];
        residualBif = @(x) bifurcation.residual(func,x,Opt,Info,Bifurcation.scaling(end));
        signDetCurrentJacobianRed = sign(detSolverJacobianRed);
        if (signDetCurrentJacobianRed*Path.signDetJRedAll(end-1)<=0) || basicFoldDetection
            bifFound=1;
            %% find exact point:
            nds = 1000;
            indBif = Path.nAll;
            bifType = NaN;
            nv = numel(Path.varAll(:,1));
            
            if signDetCurrentJacobianRed*Path.signDetJRedAll(end-1)<=0
                dss = linspace(Path.sAll(end-1)-Path.sAll(end),0,nds);
                detLeft = det(Path.getJacobianByName('previous',1:nv,1:nv));
                detRight = det(Path.getJacobianByName('last',1:nv,1:nv));
                startDsp = detRight/(detRight-detLeft)*(dss(1)-dss(end));     % estimated bifurcation point
                indStart = find(dss>startDsp,1);                                % index of estimated bifurcation point
            else
                dss = linspace(Path.sAll(end-2)-Path.sAll(end),0,nds);
                dll = interp1(Path.sAll-Path.sAll(end),abs(Path.lAll-Path.lAll(end-1)),dss,'spline');
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
                [varBifPredictor,lBifPredictor] = continuation.predictor(Path,dsp,lastJacobianRed,func,resCorr,Solver,Opt);
                dscale = aux.getDscale(Opt,struct('varAll',varBifPredictor,'lAll',lBifPredictor));
                [xBif,funBif,bifSolverExitflag,bifSolverOutput,bifSolverJacobian] = bifSolver(residualBif,[varBifPredictor;lBifPredictor],dscale);
                if bifSolverExitflag>0
                    Path.sAll = [Path.sAll(1:end-1),Path.sAll(end-1)+[norm(xBif-[Path.varAll(:,end-1);Path.lAll(end-1)]),norm(xBif-[Path.varAll(:,end-1);Path.lAll(end-1)])+norm([Path.varAll(:,end);Path.lAll(end)]-xBif)]];
                    Path.varAll = [Path.varAll(:,1:end-1),xBif(1:end-1),Path.varAll(:,end)];
                    Path.lAll = [Path.lAll(1:end-1),xBif(end),Path.lAll(end)];
                    if OptIsSet.bifAdditionalTestfunction
                        Path.bifTestValue=[Path.bifTestValue Opt.bifAdditionalTestfunction(func,xBif,Jacobian,Path,Info)];
                    end
                    % 
                    % get type of bifurcation
                    bifType = (sign(det(Path.getJacobianByName('previous')))==sign(det(Path.getJacobianByName('last')))); % 1: fold bif.; 0: branch point bif; NaN: unknown
                    %
                    % if bifType = 0 (branch point) calculate directions
                    % of paths by null() and save to Bifurcation.dirs cell array
                    if bifType == 0
                        nv = numel(Path.varAll(:,1));
                        solverJacobianRed = bifSolverJacobian(1:nv,1:nv);
                        solverJacobianLam = bifSolverJacobian(1:nv,nv+1);
                        jacRedJacLam = [solverJacobianRed, solverJacobianLam];
                        %
                        nullMatrix = null(jacRedJacLam);
                        for kDirs = 1:width(nullMatrix)
                            bifDir = nullMatrix(:,kDirs);
                            if height(Bifurcation.dirs) == 1 && isempty(Bifurcation.dirs{1,2})
                                Bifurcation.dirs(1,:) = [{indBif}, {bifDir}];
                            else
                                Bifurcation.dirs(end+1,:) = [{indBif}, {bifDir}];
                            end
                        end
                    end
                    break;
                end
            end
            Bifurcation.bif = [Bifurcation.bif,[indBif;bifType]];
            Bifurcation.flag = 1;
        end
    else
        %% error
        error('No such bifurcation-mode');
    end
    if OptIsSet.bifAdditionalTestfunction && bifFound==0
        if Opt.bifurcation.mark
            %% mark bifurcations:
            if Path.bifTestValue(end-1)*Path.bifTestValue(end)<=0
                bifFound=1;
                bifType = 2; % 2: additional; 1: fold bif.; 0: branch point bif; NaN: unknown
                Bifurcation.bif = [Bifurcation.bif,[Path.nAll;bifType]];
                Bifurcation.flag = 1;
                Bifurcation.scaling = [Bifurcation.scaling,1];
            end
        elseif Opt.bifurcation.determine || Opt.bifurcation.trace || Opt.bifurcation.parameterTrace
            %% determine bifurcation-points:
            OptBif = Opt;
            OptBif.jacobian = false;
            [bifSolver,defaultBifSolverOutput] = continuation.solver(OptBif,0);
            detSolverJacobianRed = det(lastJacobianRed);
            Bifurcation.scaling = [Bifurcation.scaling,1/detSolverJacobianRed];
            if Path.bifTestValue(end-1)*Path.bifTestValue(end)<=0
                bifFound=1;
                residualBif = @(x) bifurcation.residualAdditionalTestfunction(func,x,Opt,Jacobian,Path,Info);
                %% find exact point:
                nds = 11;
                dss = linspace(Path.sAll(end-1)-Path.sAll(end),0,nds);
                dss = dss([6,7,5,8,4,9,3,10,2,11,1]);
                indBif = Path.nAll;
                bifType = NaN;
                for ii=1:nds
                    dsp = dss(ii);
                    [varBifPredictor,lBifPredictor] = continuation.predictor(Path,dsp,lastJacobianRed,func,resCorr,Solver,Opt);
                    dscale = aux.getDscale(Opt,struct('varAll',varBifPredictor,'lAll',lBifPredictor));
                    [xBif,funBif,bifSolverExitflag,bifSolverOutput,bifSolverJacobian] = bifSolver(residualBif,[varBifPredictor;lBifPredictor],dscale);
                    if bifSolverExitflag>0
                        xBif=real(xBif);            % very small imaginary parts (probably caused by numeric differences)
                                                    % lead to errors in following code
                        Path.sAll = [Path.sAll(1:end-1),Path.sAll(end-1)+[norm(xBif-[Path.varAll(:,end-1);Path.lAll(end-1)]),norm(xBif-[Path.varAll(:,end-1);Path.lAll(end-1)])+norm([Path.varAll(:,end);Path.lAll(end)]-xBif)]];
                        Path.varAll = [Path.varAll(:,1:end-1),xBif(1:end-1),Path.varAll(:,end)];
                        Path.lAll = [Path.lAll(1:end-1),xBif(end),Path.lAll(end)];
                        % jacobian for xBif
                        Jacobianbif.solver=aux.numericJacobian(@(x)func(x(1:end-1),x(end)), xBif,'diffquot', Opt.diffquot); 
                        Path.bifTestValue=[Path.bifTestValue(1:end-1) Opt.bifAdditionalTestfunction(func,xBif,Jacobianbif,Path,Info),Path.bifTestValue(end)];
                        % 
                        % get type of bifurcation
                        bifType = 2; % 2: additional; 1: fold bif.; 0: branch point bif; NaN: unknown
                        break;
                    end
                end
                Bifurcation.bif = [Bifurcation.bif,[indBif;bifType]];
                Bifurcation.flag = 1;
            end
        else
            %% error
            error('No such bifurcation-mode');
        end
    end
end