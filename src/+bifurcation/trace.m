%% path continuation - bifurcation.trace
%
%   Institute of Dynamics and Vibration Research
%   Leibniz University Hannover
%   28.10.2020 - Alwin Förster
%   02.07.2021 - Tido Kubatschek
%
function [Path,Bifurcation] = trace(Opt,OptIsSet,Path,Bifurcation,Solver,Info,func,resCorr)
    bifTrace = Bifurcation.bif(:,Bifurcation.bif(2,:) == 0);
    nbif = numel(bifTrace(1,:));
    OptSphere = Opt;
    OptSphere = aux.seton(OptSphere,'corrector','sphere');
    OptTrace = Opt;
    OptIsSetTrace = OptIsSet;
    OptTrace.stopOnCrossing = true;
    OptIsSetTrace.stopOnCrossing = true;
    OptTrace = aux.seton(OptTrace,'bifurcation','mark');
    OptIsSetTrace.bifurcation = true;
    for i=1:nbif
        indBif = bifTrace(1,i);
        xdirsOld = [];
        xdirsTrace = [];
        dsBif = mean(diff(Path.sAll(bifTrace(1,i)+(-1:1))));
        x0 = [Path.varAll(:,bifTrace(1,i));Path.lAll(bifTrace(1,i))];
        residualBifSphere = @(x) aux.mergeResiduals(func,continuation.corrector(func,OptSphere),x,x0,dsBif,[],OptSphere);
        %% find directions of known path
        for j=1:2
            PathTrace = Path;
            if j==1
                PathTrace.varAll = Path.varAll(:,1:bifTrace(1,i));
                PathTrace.lAll = Path.lAll(1:bifTrace(1,i));
                PathTrace.sAll = Path.sAll(1:bifTrace(1,i));
            else
                PathTrace.varAll = Path.varAll(:,end:-1:bifTrace(1,i));
                PathTrace.lAll = Path.lAll(end:-1:bifTrace(1,i));
                PathTrace.sAll = abs(Path.sAll(end:-1:bifTrace(1,i))-Path.sAll(end));
            end
            [varBifPredictor,lBifPredictor] = continuation.predictor(PathTrace,dsBif,[],func,resCorr,Solver,OptSphere);
            xBifPredictor = [varBifPredictor;lBifPredictor];
            dscale = aux.getDscale(Opt,struct('varAll',varBifPredictor,'lAll',lBifPredictor));
            [xBifIj,~,solverBifExitflag] = Solver.main(residualBifSphere,xBifPredictor,dscale);
            if solverBifExitflag>0
                xdirsOld = [xdirsOld,xBifIj-x0];
                residualBifSphere = @(x) aux.deflation(residualBifSphere,xBifIj,x,OptSphere);
            end
        end
        %% unknown paths
        if ~isempty(xdirsOld)
            %% find directions of unknown paths by using random direction
            if Opt.bifRandDir
                for j=1:Opt.nBifSearch
                    dxBifPredictor = randn(numel(x0),1);
                    for ki = 1:2
                        xBifPredictor = x0+(-1)^ki*dsBif*dxBifPredictor/norm(dxBifPredictor);
                        dscale = aux.getDscale(Opt,struct('varAll',xBifPredictor(1:end-1,:),'lAll',xBifPredictor(end,:)));
                        [xBifIj,~,solverBifExitflag] = Solver.main(residualBifSphere,xBifPredictor,dscale);
                        if solverBifExitflag>0 && norm(xBifIj-x0)>=dsBif*0.99 && norm(xBifIj-x0)<=dsBif*1.01
                            xdirsTrace = [xdirsTrace,xBifIj-x0];
                            residualBifSphere = @(x) aux.deflation(residualBifSphere,xBifIj,x,OptSphere);
                        end
                    end              
                end
            end
            %% use tangent vectors - experimental!
            for j=1:(numel(Bifurcation.dirs)/2)
                [bifNum, bifDir] = Bifurcation.dirs{j,:};
                if bifNum == indBif
                    if ~isempty(bifDir)
                        for ki = 1:2
                            xBifPredictor = x0+(-1)^ki*dsBif*bifDir/norm(bifDir);
                            dscale = aux.getDscale(Opt,struct('varAll',xBifPredictor(1:end-1,:),'lAll',xBifPredictor(end,:)));
                            [xBifIj,~,solverBifExitflag] = Solver.main(residualBifSphere,xBifPredictor,dscale);
                            if solverBifExitflag>0 && norm(xBifIj-x0)>=dsBif*0.99 && norm(xBifIj-x0)<=dsBif*1.01
                                xdirsTrace = [xdirsTrace,xBifIj-x0];
                                residualBifSphere = @(x) aux.deflation(residualBifSphere,xBifIj,x,OptSphere);
                            end
                        end
                    end
                end
            end
            %% trace unknown paths
            ntrace = numel(xdirsTrace)/numel(x0);
            for j=1:ntrace
                OptTrace.direction = xdirsTrace(:,j)/norm(xdirsTrace(:,j));
                OptIsSetTrace.direction = true;
                OptTrace.l0 = x0(end);
                OptIsSetTrace.l0 = true;
                OptTrace.initialDeflationPoints = x0+xdirsOld;
                OptIsSetTrace.initialDeflationPoints = true;
                OptTrace.OptIsSet = OptIsSetTrace; % (set OptIsSet with Opt)
                [varJ,lJ,exitflagJ,BifurcationJ,sJ] = continuation(func,x0(1:end-1),Info.lStart,Info.lEnd,dsBif,'opt',OptTrace);
                if exitflagJ>0
                    Path.varAll = [Path.varAll,NaN(numel(x0)-1,1),varJ];
                    Path.lAll = [Path.lAll,NaN,lJ];
                    Path.sAll = [Path.sAll,NaN,sJ];
                    Bifurcation.bif = [Bifurcation.bif,BifurcationJ.bif];
                end
            end
        end
    end
end