%% path continuation - bifurcation.trace
%
%   Institute of Dynamics and Vibration Research
%   Leibniz University Hannover
%   28.10.2020 - Alwin Förster
%   02.07.2021 - Tido Kubatschek
%
function trace(oih,func,resCorr)
    bifTrace = oih.bifrucation.bif(:,oih.bifrucation.bif(2,:) == 0);
    nbif = numel(bifTrace(1,:));
    oihSphere = oih.copy();
    oihSphere.opt = aux.seton(oihSphere.opt,'corrector','sphere');
    oihTrace = oih.copy();
    oihTrace.opt.stopOnCrossing = true;
    oihTrace.optIsSet.stopOnCrossing = true;
    oihTrace.opt = aux.seton(oihTrace.opt,'bifurcation','mark');
    oihTrace.optIsSet.bifurcation = true;
    for i=1:nbif
        indBif = bifTrace(1,i);
        xdirsOld = [];
        xdirsTrace = [];
        dsBif = mean(diff(oih.path.sAll(bifTrace(1,i)+(-1:1))));
        x0 = [oih.path.varAll(:,bifTrace(1,i));oih.path.lAll(bifTrace(1,i))];
        residualBifSphere = @(x) aux.mergeResiduals(func,continuation.corrector(func,oihSphere),x,x0,dsBif,[],oihSphere);
        %% find directions of known path
        for j=1:2
            if j==1
                varAllTemp = oih.path.varAll(:,1:bifTrace(1,i));
                lAllTemp = oih.path.lAll(1:bifTrace(1,i));
            else
                varAllTemp = oih.path.varAll(:,end:-1:bifTrace(1,i));
                lAllTemp = oih.path.lAll(end:-1:bifTrace(1,i));
            end
            [varBifPredictor,lBifPredictor] = continuation.predictor(oihSphere,dsBif,[],func,resCorr);
            xBifPredictor = [varBifPredictor;lBifPredictor];
            dscale = aux.getDscale(oih,struct('varAll',varBifPredictor,'lAll',lBifPredictor));
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
                for j=1:oih.opt.nBifSearch
                    dxBifPredictor = randn(numel(x0),1);
                    for ki = 1:2
                        xBifPredictor = x0+(-1)^ki*dsBif*dxBifPredictor/norm(dxBifPredictor);
                        dscale = aux.getDscale(oih,struct('varAll',xBifPredictor(1:end-1,:),'lAll',xBifPredictor(end,:)));
                        [xBifIj,~,solverBifExitflag] = oih.solver.main(residualBifSphere,xBifPredictor,dscale);
                        if solverBifExitflag>0 && norm(xBifIj-x0)>=dsBif*0.99 && norm(xBifIj-x0)<=dsBif*1.01
                            xdirsTrace = [xdirsTrace,xBifIj-x0];
                            residualBifSphere = @(x) aux.deflation(residualBifSphere,xBifIj,x,oih.opt.jacobian);
                        end
                    end              
                end
            end
            %% use tangent vectors - experimental!
            for j=1:(numel(oih.bifrucation.dirs)/2)
                [bifNum, bifDir] = oih.bifrucation.dirs{j,:};
                if bifNum == indBif
                    if ~isempty(bifDir)
                        for ki = 1:2
                            xBifPredictor = x0+(-1)^ki*dsBif*bifDir/norm(bifDir);
                            dscale = aux.getDscale(oih,struct('varAll',xBifPredictor(1:end-1,:),'lAll',xBifPredictor(end,:)));
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
            for j=1:ntrace
                oihTrace.opt.direction = xdirsTrace(:,j)/norm(xdirsTrace(:,j));
                oihTrace.optIsSet.direction = true;
                oihTrace.opt.l0 = x0(end);
                oihTrace.optIsSet.l0 = true;
                oihTrace.opt.initialDeflationPoints = x0+xdirsOld;
                oihTrace.optIsSet.initialDeflationPoints = true;
                [varJ,lJ,exitflagJ,bifStructJ,sJ] = continuation(func,x0(1:end-1),oih.info.lStart,oih.info.lEnd,dsBif,'opt',oihTrace.opt);
                if exitflagJ>=0
                    varAllTemp = [oih.path.varAll,NaN(numel(x0)-1,1),varJ];
                    lAllTemp = [oih.path.lAll,NaN,lJ];
                    bifTemp = [oih.bifrucation.bif,bifStructJ.bif];
                    oih.path.overwrite(varAllTemp,lAllTemp,bifTemp,true);
                end
            end
        end
    end
end