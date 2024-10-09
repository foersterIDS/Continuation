%% path continuation - bifurcation.trace
%
%   Institute of Dynamics and Vibration Research
%   Leibniz University Hannover
%   28.10.2020 - Alwin Förster
%   02.07.2021 - Tido Kubatschek
%
function trace(oih,func,resCorr)
    bifTrace = oih.bifurcation.bif(:,oih.bifurcation.bif(2,:) == 0);
    nbif = numel(bifTrace(1,:));
    oihSphere = oih.copy();
    oihSphere.opt = aux.seton(oihSphere.opt,'corrector','sphere');
    oihTrace = oih.copy();
    oihTrace.opt.stopOnCrossing = true;
    oihTrace.optIsSet.stopOnCrossing = true;
    oihTrace.opt = aux.seton(oihTrace.opt,'bifurcation','mark');
    oihTrace.optIsSet.bifurcation = true;
    for ii=1:nbif
        indBif = bifTrace(1,ii);
        xdirsOld = [];
        xdirsTrace = [];
        dsBif = mean(diff(oih.path.sAll(bifTrace(1,ii)+(-1:1))));
        x0 = [oih.path.varAll(:,bifTrace(1,ii));oih.path.lAll(bifTrace(1,ii))];
        residualBifSphere = @(x) aux.mergeResiduals(func,continuation.corrector(func,oihSphere),x,x0,dsBif,[],oihSphere);
        %% find directions of known path
        for jj=1:2
            if jj==1
                varAllTemp = oih.path.varAll(:,1:bifTrace(1,ii));
                lAllTemp = oih.path.lAll(1:bifTrace(1,ii));
            else
                varAllTemp = oih.path.varAll(:,end:-1:bifTrace(1,ii));
                lAllTemp = oih.path.lAll(end:-1:bifTrace(1,ii));
            end
            [varBifPredictor,lBifPredictor] = continuation.predictor(oihSphere,dsBif,[],func,resCorr);
            xBifPredictor = [varBifPredictor;lBifPredictor];
            dscale = aux.getDscale(oih,xBifPredictor);
            [xBifIj,~,solverBifExitflag] = oih.solver.main(residualBifSphere,xBifPredictor,dscale);
            if solverBifExitflag>0
                xdirsOld = [xdirsOld,xBifIj-x0];
                residualBifSphere = @(x) aux.deflation(residualBifSphere,xBifIj,x,oih.opt.jacobian);
            end
        end
        %% unknown paths
        if ~isempty(xdirsOld)
            %% find directions of unknown paths by using random direction
            if oih.opt.bifRandDir
                for jj=1:oih.opt.nBifSearch
                    dxBifPredictor = randn(numel(x0),1);
                    for ki = 1:2
                        xBifPredictor = x0+(-1)^ki*dsBif*dxBifPredictor/norm(dxBifPredictor);
                        dscale = aux.getDscale(oih,xBifPredictor);
                        [xBifIj,~,solverBifExitflag] = oih.solver.main(residualBifSphere,xBifPredictor,dscale);
                        if solverBifExitflag>0 && norm(xBifIj-x0)>=dsBif*0.99 && norm(xBifIj-x0)<=dsBif*1.01
                            xdirsTrace = [xdirsTrace,xBifIj-x0];
                            residualBifSphere = @(x) aux.deflation(residualBifSphere,xBifIj,x,oih.opt.jacobian);
                        end
                    end              
                end
            end
            %% use tangent vectors - experimental!
            for jj=1:(numel(oih.bifurcation.dirs)/2)
                [bifNum, bifDir] = oih.bifurcation.dirs{jj,:};
                if bifNum == indBif
                    if ~isempty(bifDir)
                        for ki = 1:2
                            xBifPredictor = x0+(-1)^ki*dsBif*bifDir/norm(bifDir);
                            dscale = aux.getDscale(oih,xBifPredictor);
                            [xBifIj,~,solverBifExitflag] = oih.solver.main(residualBifSphere,xBifPredictor,dscale);
                            if solverBifExitflag>0 && norm(xBifIj-x0)>=dsBif*0.99 && norm(xBifIj-x0)<=dsBif*1.01
                                xdirsTrace = [xdirsTrace,xBifIj-x0];
                                residualBifSphere = @(x) aux.deflation(residualBifSphere,xBifIj,x,oih.opt.jacobian);
                            end
                        end
                    end
                end
            end
            %% trace unknown paths
            ntrace = numel(xdirsTrace)/numel(x0);
            optTrace = oihTrace.opt;
            fnOiS = fieldnames(oihTrace.optIsSet);
            optTrace = rmfield(optTrace,fnOiS(~cell2mat(struct2cell(oihTrace.optIsSet)))); 
            for jj=1:ntrace
                optTrace.direction = xdirsTrace(:,jj)/norm(xdirsTrace(:,jj));
                optTrace.l0 = x0(end);
                optTrace.initialDeflationPoints = x0+xdirsOld;
                [varJ,lJ,exitflagJ,bifStructJ,sJ] = continuation(func,x0(1:end-1),oih.info.lStart,oih.info.lEnd,dsBif,'opt',optTrace);
                if exitflagJ>=0
                    varAllTemp = [oih.path.varAll,NaN(numel(x0)-1,1),varJ];
                    lAllTemp = [oih.path.lAll,NaN,lJ];
                    bifTemp = [oih.bifurcation.bif,bifStructJ.bif];
                    oih.path.overwrite(varAllTemp,lAllTemp,bifTemp,true);
                end
            end
        end
    end
end