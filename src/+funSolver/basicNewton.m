%% path continuation - funSolver.basicNewton
%
%   Institute of Dynamics and Vibration Research
%   Leibniz University Hannover
%   19.05.2020 - Alwin Förster
%   16.09.2020 - Tido Kubatschek 
%
function [x,fval,exitflag,output,jacobian] = basicNewton(fun,x0,dscale,outputFlag,Opt)
    maxStep = Opt.solverMaxIterations;
    if outputFlag
        global solverStepsizes;
        solverStepsizes = [0, 0];
    end
    nSteps = 0;
    exitflag = -1;
    xi = x0;
    xim1 = x0;
    xiSc = x0./dscale;
    for i=1:maxStep
        if Opt.jacobian
            [fi,Ji] = fun(xi);
        else
            fi = fun(xi);
            Ji = aux.numericJacobian(fun,xi,'centralValue',fi,'diffquot',Opt.diffquot,'diffStep',oih.opt.diffStep);
        end
        JiSc = Ji*diag(dscale);
        nf = length(fi);
        xip1Sc = xiSc-JiSc(1:nf,1:nf)\fi;
        %% end loop
        xiSc = xip1Sc;
        xi = xip1Sc.*dscale;
        nSteps = nSteps+1;
        absFi = sqrt(fi'*fi);
        absDxi = norm(xi-xim1)/norm(xi);
        if outputFlag
            solverStepsizes = [solverStepsizes; i, norm(xi-xim1)];
        end
        if absFi<Opt.solverTol
            exitflag = 2;
            break;
        elseif absDxi<Opt.solverTol
            exitflag = 1;
            break;
        end
        xim1 = xi;
    end
    %% make output
    x = xi;
    if nSteps==maxStep
        % TODO: hier muss es noch mehr geben
        exitflag = 0;
    end
    if nargout>1
        if Opt.jacobian
            [fi,Ji] = fun(xi);
        else
            fi = fun(xi);
        end
        fval = fi;
        if nargout>3
            output = struct('iterations',nSteps,...
                'algorithm','basic-newton',...
                'tolerance',absFi);
            if nargout>4
                if ~Opt.jacobian
                    Ji = aux.numericJacobian(fun,xi,'centralValue',fi,'diffquot',Opt.diffquot,'diffStep',oih.opt.diffStep);
                end
                jacobian = Ji;
            end
        end
    end
end