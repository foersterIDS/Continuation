%% path continuation - corrector.adapt
%
%   Institute of Dynamics and Vibration Research
%   Leibniz University Hannover
%   20.01.2022 - Alwin Förster
%
function [Do,Opt,corrInfo] = adapt(Do,Opt,Path,Solver,func,xPredictor,dscale,lastJacobian,ds)
    corrInfo = [];
    if aux.ison(Opt.adaptCorrector)
        if Opt.corrector.sphere
            if Solver.exitflag>0
                if Opt.adaptCorrector.solve
                    OptTemp = aux.seton(Opt,'corrector','orthogonal');
                    resCorrTemp = continuation.corrector(func,OptTemp);
                    residualTemp = @(x) aux.mergeResiduals(func,resCorrTemp,x,[Path.varAll;Path.lAll],ds,lastJacobian,OptTemp);
                    [~,~,solverExitflagTemp,solverOutputTemp] = Solver.main(residualTemp,xPredictor,dscale);
                    if solverExitflagTemp>0 && solverOutputTemp.iterations<Solver.output.iterations
                        corrInfo = 'orthogonal';
                        Do.changeCorrector = true;
                    end
                end
            else
                corrInfo = 'orthogonal';
                Do.changeCorrector = true;
            end
        elseif Opt.corrector.orthogonal
            if Solver.exitflag>0
                if Opt.adaptCorrector.solve
                    OptTemp = aux.seton(Opt,'corrector','sphere');
                    resCorrTemp = continuation.corrector(OptTemp);
                    residualTemp = @(x) aux.mergeResiduals(func,resCorrTemp,x,[Path.varAll;Path.lAll],ds,lastJacobian,OptTemp);
                    [~,~,solverExitflagTemp,solverOutputTemp] = Solver.main(residualTemp,xPredictor,dscale);
                    if solverExitflagTemp>0 && solverOutputTemp.iterations<Solver.output.iterations
                        corrInfo = 'sphere';
                        Do.changeCorrector = true;
                    end
                end
            else
                corrInfo = 'sphere';
                Do.changeCorrector = true;
            end
        else
            corrInfo = 'orthogonal';
            Do.changeCorrector = true;
        end
    end
end