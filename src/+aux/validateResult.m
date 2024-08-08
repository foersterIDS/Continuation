%% path continuation - aux.validateResult
%
%   Institute of Dynamics and Vibration Research
%   Leibniz University Hannover
%   08.05.2020 - Alwin Förster
%
function [invPoiStr,Counter,Do,Is,Opt] = validateResult(xSolution,funSolution,Path,ds,Solver,funPredictor,sPredictor,Do,Bifurcation,Info,Is,Counter,Plot,Opt)
    %% automated validation
    %
    Is.reverse = false;
    Is.catch = 0;
    invPoiStr = '                          ';
    nPath = Path.nAll;
    if Solver.exitflag>0
        try
            if ~Opt.checkResidual || (norm(funSolution)<=Opt.solverTol*10)
                xi = [Path.varAll(:,end);Path.lAll(end)];
                normXsXi = norm(xSolution-xi);
                if nPath>1 && Opt.alphaReverseAutoMode
                    %% alphaReverseAutoMode
                    dsHist = sqrt(sum(([Path.varAll;Path.lAll]-[Path.varAll(:,end);Path.lAll(end)]).^2));
                    oobHist = find(dsHist>2*ds);
                    if ~isempty(oobHist) && oobHist(end)<(nPath-1)
                        idxHist = oobHist(end):(nPath-1);
                        xHist = [Path.varAll(:,idxHist);Path.lAll(idxHist)];
                        alphaHist = acos(((xHist(:,1:(end-1))-xi)'*(xHist(:,end)-xi))./(sqrt(diag((xHist(:,1:(end-1))-xi)'*(xHist(:,1:(end-1))-xi)))*sqrt((xHist(:,end)-xi)'*(xHist(:,end)-xi))));
                        Opt.alphaReverse = min(max(pi-4*max(alphaHist),pi/64),3/4*pi);
                    end
                end
                if ((normXsXi>=Opt.dsTol(1)*ds && normXsXi<=Opt.dsTol(2)*ds || nPath==1) && normXsXi<=Opt.dsTol(2)*Opt.dsMax) || Do.convergeToTarget || Opt.corrector.unique
                    if nPath==1
                        if numel(Opt.direction)==1 && sign(xSolution(end)-Path.lAll(end))==sign(Opt.direction)
                            Is.valid = true;
                        else
                            alpha = acos(((xSolution-xi)'*Opt.direction)/(sqrt((xSolution-xi)'*(xSolution-xi))*sqrt(Opt.direction'*Opt.direction)));
                            if alpha<Opt.alphaReverse
                                Is.valid = true;
                            else
                                Is.valid = false;
                                Is.reverse = true;
                                invPoiStrTemp = sprintf('(alpha: %.2f)',alpha);
                                invPoiStr(1:numel(invPoiStrTemp)) = invPoiStrTemp;
                            end
                        end
                    else
                        xim1 = [Path.varAll(:,end-1);Path.lAll(end-1)];
                        alpha = acos(((xSolution-xi)'*(xi-xim1))/(sqrt((xSolution-xi)'*(xSolution-xi))*sqrt((xi-xim1)'*(xi-xim1))));
                        if alpha<Opt.alphaReverse
                            if aux.ison(Opt.bifurcation) && ~Opt.bifurcation.mark
                                if sign(det(Solver.jacobian(1:Info.nv,1:Info.nv)))*Path.signDetJRedAll(:,end)<0 && Counter.bifStepsizeRed<2
                                    Is.valid = false;
                                    Do.stepbackManually = true;
                                    Counter.bifStepsizeRed = Counter.bifStepsizeRed+1;
                                else
                                    Is.valid = true;
                                end
                            else
                                Is.valid = true;
                            end
                        else
                            if ~isempty(Bifurcation.bif) && Bifurcation.bif(1,end)==nPath
                                Is.valid = true;
                            else
                                Is.valid = false;
                                Is.reverse = true;
                                invPoiStrTemp = sprintf('(alpha: %.2f)',alpha);
                                invPoiStr(1:numel(invPoiStrTemp)) = invPoiStrTemp;
                            end
                        end
                    end
                else
                    Is.valid = false;
                    invPoiStr(1:25) = '(dx outside of ds-bounds)';
                end
            else
                Is.valid = false;
                invPoiStr(1:18) = '(fun=0 not solved)';
            end
        catch
            Is.valid = false;
            Is.catch = 1;
            invPoiStr(1:24) = '(error while validating)';
        end
    else
        Is.valid = false;
        invPoiStr(1:26) = '(Solver does not converge)';
    end
    %
    %% enforceDsMax
    %
    if Is.valid && Opt.enforceDsMax
        Is.valid = logical(heaviside(min(Opt.dsMax-(xSolution-xi))));
        if ~Is.valid
            invPoiStr(1:17) = '(dsMax exceeded)';
        end
    end
    %
    %% approve manually
    %
    if Is.valid && ((islogical(Opt.approveManually) && Opt.approveManually) || (~islogical(Opt.approveManually) && (xSolution(end)-Opt.approveManually)*(Path.lAll(end)-Opt.approveManually)<=0))
        Opt.approveManually = true;
        if nPath>1
            try
                [Plot, Opt] = livePlot(Opt, Info, Path, ds, ds, Solver.output.iterations, Counter, funPredictor, sPredictor, Plot, Bifurcation);
            catch
                aux.printLine(Opt,'--> The plot update for approval has failed.\n');
            end
        end
        prompt = sprintf('------> approve point at l = %.4e (y/n): ',Path.lAll(end));
        correctInput = false;
        while ~correctInput
            inputString = input(prompt,'s');
            if strcmp(inputString,'y')
                correctInput = true;
            elseif strcmp(inputString,'n')
                correctInput = true;
                Is.valid = false;
                invPoiStr(1:19) = '(rejected manually)';
            elseif strcmp(inputString,'off')
                correctInput = true;
                Opt.approveManually = false;
            elseif strcmp(inputString,'exit')
                correctInput = true;
                Is.valid = false;
                Do.stopManually = true;
                invPoiStr(1:19) = '(rejected manually)';
            elseif ~isnan(str2double(inputString))
                correctInput = true;
                Opt.approveManually = str2double(inputString);
            else
                aux.printLine(Opt,'------> Enter ''y'' for yes or ''n'' for no! (Deactivate with ''off'', set double limit or leave using ''exit'')\n');
            end
        end
    end
    %
    %% finish
    %
    if Is.valid
        if Do.stepbackManually
            Do.stepbackManually = false;
        else
            Counter.bifStepsizeRed = 0;
        end
    end
    %
end