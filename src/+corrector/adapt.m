%% path continuation - corrector.adapt
%
%   Institute of Dynamics and Vibration Research
%   Leibniz University Hannover
%   20.01.2022 - Alwin Förster
%
function [corrInfo] = adapt(oih,func,xPredictor,dscale,lastJacobian,ds)
    corrInfo = [];
    if aux.ison(oih.opt.adaptCorrector)
        if oih.opt.corrector.sphere
            if oih.solver.exitflag>0
                if oih.opt.adaptCorrector.solve
                    [~,correctorTemp] = ison(oih.opt.corrector);
                    oih.opt.corrector = aux.seton(ohi.opt,'corrector','orthogonal');
                    resCorrTemp = continuation.corrector(func,oih);
                    residualTemp = @(x) aux.mergeResiduals(func,resCorrTemp,x,oih.path.xAll,ds,lastJacobian,oih);
                    [~,~,solverExitflagTemp,solverOutputTemp] = oih.solver.main(residualTemp,xPredictor,dscale);
                    oih.opt.corrector = aux.seton(ohi.opt,'corrector',correctorTemp);
                    if solverExitflagTemp>0 && solverOutputTemp.iterations<oih.solver.output.iterations
                        corrInfo = 'orthogonal';
                        oih.do.changeCorrector = true;
                    end
                end
            else
                corrInfo = 'orthogonal';
                oih.do.changeCorrector = true;
            end
        elseif oih.opt.corrector.orthogonal
            if oih.solver.exitflag>0
                if oih.opt.adaptCorrector.solve
                    [~,correctorTemp] = ison(oih.opt.corrector);
                    oih.opt.corrector = aux.seton(oih.opt,'corrector','sphere');
                    resCorrTemp = continuation.corrector(func,oih);
                    residualTemp = @(x) aux.mergeResiduals(func,resCorrTemp,x,oih.path.xAll,ds,lastJacobian,OptTemp);
                    [~,~,solverExitflagTemp,solverOutputTemp] = oih.solver.main(residualTemp,xPredictor,dscale);
                    oih.opt.corrector = aux.seton(ohi.opt,'corrector',correctorTemp);
                    if solverExitflagTemp>0 && solverOutputTemp.iterations<oih.solver.output.iterations
                        corrInfo = 'sphere';
                        oih.do.changeCorrector = true;
                    end
                end
            else
                corrInfo = 'sphere';
                oih.do.changeCorrector = true;
            end
        else
            corrInfo = 'orthogonal';
            oih.do.changeCorrector = true;
        end
    end
end