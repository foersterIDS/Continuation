%% path continuation - continuation.solver
%
%   Institute of Dynamics and Vibration Research
%   Leibniz University Hannover
%   08.05.2020 - Alwin Förster
%
function [solver,predictorSolver,numJacSolver,defaultSolverOutput] = solver(Opt,outputFlag)
    predictorSolver = [];
    if outputFlag
        outfun = @outputFun;
    else
        outfun = @doNothing;
    end
    if Opt.solver.fsolve
        %% fsolve
        if Opt.jacobian
            options = optimoptions('fsolve','display','off','SpecifyObjectiveGradient',true,'FunctionTolerance',Opt.solverTol,'StepTolerance',Opt.solverTol,'OptimalityTolerance',Opt.solverTol,'MaxIterations',Opt.solverMaxIterations,'OutputFcn', outfun);
        else
            if Opt.diffquot.central
                options = optimoptions('fsolve','display','off','SpecifyObjectiveGradient',false,'FiniteDifferenceType','central','FunctionTolerance',Opt.solverTol,'StepTolerance',Opt.solverTol,'OptimalityTolerance',Opt.solverTol,'MaxIterations',Opt.solverMaxIterations,'OutputFcn', outfun);
            else
                options = optimoptions('fsolve','display','off','SpecifyObjectiveGradient',false,'FiniteDifferenceType','forward','FunctionTolerance',Opt.solverTol,'StepTolerance',Opt.solverTol,'OptimalityTolerance',Opt.solverTol,'MaxIterations',Opt.solverMaxIterations,'OutputFcn', outfun);
            end
        end
        optfun = @(dscale) optimoptions(options,'TypicalX',dscale);
        solver = @(fun,x0,dscale) fsolve(fun,x0,optfun(dscale));
        if Opt.predictorSolver
            predictorOptions = optimoptions('fsolve','display','off','SpecifyObjectiveGradient',true,'FunctionTolerance',Opt.solverTol,'StepTolerance',Opt.solverTol,'OptimalityTolerance',Opt.solverTol,'MaxIterations',Opt.solverMaxIterations);
            predictorSolver = @(funPredictor,s0) fsolve(funPredictor,s0,predictorOptions);
        end
        optNumJac = optimoptions(options,'SpecifyObjectiveGradient',false);
        numJacSolver = @(fun,x0,dscale) fsolve(fun,x0,optNumJac);
    elseif Opt.solver.lsqnonlin
        %% lsqnonlin
        if Opt.jacobian
            options = optimoptions('lsqnonlin','display','off','SpecifyObjectiveGradient',true,'FunctionTolerance',Opt.solverTol,'MaxIterations',Opt.solverMaxIterations,'OutputFcn', outfun);
        else
            if Opt.diffquot.central
                options = optimoptions('lsqnonlin','display','off','SpecifyObjectiveGradient',false,'FiniteDifferenceType','central','FunctionTolerance',Opt.solverTol,'MaxIterations',Opt.solverMaxIterations,'OutputFcn', outfun);
            else
                options = optimoptions('lsqnonlin','display','off','SpecifyObjectiveGradient',false,'FiniteDifferenceType','forward','FunctionTolerance',Opt.solverTol,'MaxIterations',Opt.solverMaxIterations,'OutputFcn', outfun);
            end
        end
        optfun = @(dscale) optimoptions(options,'TypicalX',dscale);
        solver = @(fun,x0,dscale) funSolver.sLsqnonlin(fun,x0,optfun(dscale));
        if Opt.predictorSolver
            predictorOptions = optimoptions('lsqnonlin','display','off','SpecifyObjectiveGradient',true,'FunctionTolerance',Opt.solverTol,'MaxIterations',Opt.solverMaxIterations);
            predictorSolver = @(funPredictor,s0) funSolver.sLsqnonline(funPredictor,s0,predictorOptions);
        end
        optNumJac = optimoptions(options,'SpecifyObjectiveGradient',false);
        numJacSolver = @(fun,x0,dscale) funSolver.sLsqnonline(fun,x0,optNumJac);
    elseif Opt.solver.newton
        %% basic newton solver
        solver = @(fun,x0,dscale) funSolver.basicNewton(fun,x0,dscale,outputFlag,Opt);
        if Opt.predictorSolver
            predictorOpt = Opt;
            predictorOpt.jacobian = true;
            predictorSolver = @(funPredictor,s0) funSolver.basicNewton(funPredictor,s0,1,predictorOpt);
        end
        optNumJac = Opt;
        optNumJac.jacobian = false;
        numJacSolver = @(fun,x0,dscale) funSolver.basicNewton(fun,x0,1,optNumJac);
    else
        %% error
        error('No such solver');
    end
    defaultSolverOutput = struct('iterations',inf,...
                                   'algorithm','non',...
                                   'tolerance',inf);
end

function stop = outputFun(x,optimValues,state)
    stop = false;
    global solverStepsizes;
    switch state
        case 'init'
            solverStepsizes = [];
        case 'iter'
            iter = optimValues.iteration;               % Iteration
            normd = optimValues.stepsize;               % Norm of step
            solverStepsizes = [solverStepsizes; iter normd];
    end
end

function stop = doNothing(x,optimValues,state)
    stop = false;
end