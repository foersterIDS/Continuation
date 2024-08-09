%% path continuation - aux.validateResult
%
%   Institute of Dynamics and Vibration Research
%   Leibniz University Hannover
%   08.05.2020 - Alwin Förster
%
function [invPoiStr] = validateResult(xSolution,funSolution,oih,ds,funPredictor,sPredictor)
    %% automated validation
    %
    oih.is.reverse = false;
    oih.is.catch = 0;
    invPoiStr = '                          ';
    nPath = oih.path.nAll;
    if oih.solver.exitflag>0
        try
            if ~oih.opt.checkResidual || (norm(funSolution)<=oih.opt.solverTol*10)
                xi = [oih.path.varAll(:,end);oih.path.lAll(end)];
                normXsXi = norm(xSolution-xi);
                if nPath>1 && oih.opt.alphaReverseAutoMode
                    %% alphaReverseAutoMode
                    dsHist = sqrt(sum(([oih.path.varAll;oih.path.lAll]-[oih.path.varAll(:,end);oih.path.lAll(end)]).^2));
                    oobHist = find(dsHist>2*ds);
                    if ~isempty(oobHist) && oobHist(end)<(nPath-1)
                        idxHist = oobHist(end):(nPath-1);
                        xHist = [oih.path.varAll(:,idxHist);oih.path.lAll(idxHist)];
                        alphaHist = acos(((xHist(:,1:(end-1))-xi)'*(xHist(:,end)-xi))./(sqrt(diag((xHist(:,1:(end-1))-xi)'*(xHist(:,1:(end-1))-xi)))*sqrt((xHist(:,end)-xi)'*(xHist(:,end)-xi))));
                        oih.opt.alphaReverse = min(max(pi-4*max(alphaHist),pi/64),3/4*pi);
                    end
                end
                if ((normXsXi>=oih.opt.dsTol(1)*ds && normXsXi<=oih.opt.dsTol(2)*ds || nPath==1) && normXsXi<=oih.opt.dsTol(2)*oih.opt.dsMax) || oih.do.convergeToTarget || oih.opt.corrector.unique
                    if nPath==1
                        if numel(oih.opt.direction)==1 && sign(xSolution(end)-oih.path.lAll(end))==sign(oih.opt.direction)
                            oih.is.valid = true;
                        else
                            alpha = acos(((xSolution-xi)'*oih.opt.direction)/(sqrt((xSolution-xi)'*(xSolution-xi))*sqrt(oih.opt.direction'*oih.opt.direction)));
                            if alpha<oih.opt.alphaReverse
                                oih.is.valid = true;
                            else
                                oih.is.valid = false;
                                oih.is.reverse = true;
                                invPoiStrTemp = sprintf('(alpha: %.2f)',alpha);
                                invPoiStr(1:numel(invPoiStrTemp)) = invPoiStrTemp;
                            end
                        end
                    else
                        xim1 = [oih.path.varAll(:,end-1);oih.path.lAll(end-1)];
                        alpha = acos(((xSolution-xi)'*(xi-xim1))/(sqrt((xSolution-xi)'*(xSolution-xi))*sqrt((xi-xim1)'*(xi-xim1))));
                        if alpha<oih.opt.alphaReverse
                            if aux.ison(oih.opt.bifurcation) && ~oih.opt.bifurcation.mark
                                if sign(det(oih.solver.jacobian(1:oih.info.nv,1:oih.info.nv)))*oih.path.signDetJRedAll(:,end)<0 && oih.counter.bifStepsizeRed<2
                                    oih.is.valid = false;
                                    oih.do.stepbackManually = true;
                                    oih.counter.bifStepsizeRed = oih.counter.bifStepsizeRed+1;
                                else
                                    oih.is.valid = true;
                                end
                            else
                                oih.is.valid = true;
                            end
                        else
                            if ~isempty(oih.bifurcation.bif) && oih.bifurcation.bif(1,end)==nPath
                                oih.is.valid = true;
                            else
                                oih.is.valid = false;
                                oih.is.reverse = true;
                                invPoiStrTemp = sprintf('(alpha: %.2f)',alpha);
                                invPoiStr(1:numel(invPoiStrTemp)) = invPoiStrTemp;
                            end
                        end
                    end
                else
                    oih.is.valid = false;
                    invPoiStr(1:25) = '(dx outside of ds-bounds)';
                end
            else
                oih.is.valid = false;
                invPoiStr(1:18) = '(fun=0 not solved)';
            end
        catch
            oih.is.valid = false;
            oih.is.catch = 1;
            invPoiStr(1:24) = '(error while validating)';
        end
    else
        oih.is.valid = false;
        invPoiStr(1:26) = '(Solver does not converge)';
    end
    %
    %% enforceDsMax
    %
    if oih.is.valid && oih.opt.enforceDsMax
        oih.is.valid = logical(heaviside(min(oih.opt.dsMax-(xSolution-xi))));
        if ~oih.is.valid
            invPoiStr(1:17) = '(dsMax exceeded)';
        end
    end
    %
    %% approve manually
    %
    if oih.is.valid && ((islogical(oih.opt.approveManually) && oih.opt.approveManually) || (~islogical(oih.opt.approveManually) && (xSolution(end)-oih.opt.approveManually)*(oih.path.lAll(end)-oih.opt.approveManually)<=0))
        oih.opt.approveManually = true;
        if nPath>1
            try
                livePlot(oih, ds, ds, oih.solver.output.iterations, funPredictor, sPredictor);
            catch
                aux.printLine(oih,'--> The plot update for approval has failed.\n');
            end
        end
        prompt = sprintf('------> approve point at l = %.4e (y/n): ',oih.path.lAll(end));
        correctInput = false;
        while ~correctInput
            inputString = input(prompt,'s');
            if strcmp(inputString,'y')
                correctInput = true;
            elseif strcmp(inputString,'n')
                correctInput = true;
                oih.is.valid = false;
                invPoiStr(1:19) = '(rejected manually)';
            elseif strcmp(inputString,'off')
                correctInput = true;
                oih.opt.approveManually = false;
            elseif strcmp(inputString,'exit')
                correctInput = true;
                oih.is.valid = false;
                oih.do.stopManually = true;
                invPoiStr(1:19) = '(rejected manually)';
            elseif ~isnan(str2double(inputString))
                correctInput = true;
                oih.opt.approveManually = str2double(inputString);
            else
                aux.printLine(oih,'------> Enter ''y'' for yes or ''n'' for no! (Deactivate with ''off'', set double limit or leave using ''exit'')\n');
            end
        end
    end
    %
    %% finish
    %
    if oih.is.valid
        if oih.do.stepbackManually
            oih.do.stepbackManually = false;
        else
            oih.counter.bifStepsizeRed = 0;
        end
    end
    %
end