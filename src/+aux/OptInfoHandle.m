%% path continuation - aux.OptInfoHandle (class)
%
%   Institute of Dynamics and Vibration Research
%   Leibniz University Hannover
%   09.08.2024 - Alwin Förster
%
classdef OptInfoHandle < handle
    %OptInfoHandle
    %   Handle class that contains structs and objects with options and
    %   information.

    properties
        bifurcation
        counter
        do
        info
        infoOut
        initial
        is
        opt
        optIsSet
        path
        plot
        remove
        solver
        stepsizeOptions
    end

    methods
        %% constructor
        function obj = OptInfoHandle(opt,optIsSet)
            obj.opt = opt;
            obj.optIsSet = optIsSet;
        end

        %% general functions
        function oihOut = copy(obj)
            oihOut = aux.OptInfoHandle(obj.opt,obj.optIsSet);
            props = properties(obj);
            for ii=1:numel(props)
                if isstruct(obj.(props{ii}))
                    oihOut.(props{ii}) = obj.(props{ii});
                else
                    oihOut.(props{ii}) = obj.(props{ii}).copy();
                end
            end
        end

        function initializeStructsAndClasses(obj,var0,lStart,lEnd,ds0)
            %% Bifurcation
            %
            obj.bifurcation = struct('bif',zeros(2,0),...
                                     'dirs',{cell(1,2)},...
                                     'flag',0,...
                                     'scaling',[]);
            %
            %% Counter
            %
            obj.counter = struct('bifStepsizeRed',0,...
                                 'catch',0,...
                                 'catchOld',0,...
                                 'closedCurve',0,...
                                 'error',0,...
                                 'event',0,...
                                 'loop',0,...
                                 'remove',0,...
                                 'step',0,...
                                 'validStepback',0);
            %
            %% Do
            %
            obj.do = struct('changeCorrector',false,...
                            'continuation',false,...
                            'convergeToTarget',false,...
                            'convergeToxTarget',false,...
                            'deflate',false,...
                            'homotopy',false,...
                            'loop',false,...
                            'remove',false,...
                            'stepback',false,...
                            'stepbackManually',false,...
                            'stopManually',false,...
                            'suspend',false);
        	%
            %% Info
            %
            obj.info = struct('biDirRuns',0,...
                              'validJacobian',true,...
                              'ds0',ds0,...
                              'exitflag',-1,...
                              'exitMsg','',...
                              'nv',numel(var0),...
                              'nl',1,...
                              'lStart',lStart,...
                              'lEnd',lEnd,...
            				  't0',tic,...
                              'var0',var0,...
                              'finalSolutionPoint',false);
        	%
            %% InfoOut
            %
            obj.infoOut = struct('numberOfSteps',0,...
                                 'numberOfInvalidPoints',0);
            %
            %% Initial
            %
            obj.initial = struct('dsMax',obj.opt.dsMax,...
                                 'dsMin',obj.opt.dsMin,...
                                 'lStart',lStart,...
                                 'lEnd',lEnd);
            %
            %% Is
            %
            obj.is = struct('catch',false,...
                            'currentJacobian',false,...
                            'reverse',false,...
                            'valid',false);
            %
            %% Path
            %
            if obj.optIsSet.lMult0
                obj.info.nl = numel(obj.opt.lMult0);
                obj.opt.direction = 1;
            else
                obj.info.nl = 1;
            end
            obj.path = continuation.Path(obj.info.nv,obj.info.nl,obj);
        	%
            %% Plot
            %
            if (aux.ison(obj.opt.plot) && ~obj.opt.plot.detail) || ~aux.ison(obj.opt.plot)
                obj.plot = struct('fig',[],...
                                  'pl',[],...
                                  'plCurr',[]);
            else
                obj.plot = struct('fig',[],...
                                  'pl',[],...
                                  'plCor',[],...
                                  'plCorAssist',[],...
                                  'plCurr',[],...
                                  'plDet',[],...
                                  'plIt',[],...
                                  'plPre',[],...
                                  'plS',[]);
            end
        	%
            %% Remove
            %
            obj.remove = struct('ds',NaN,...
                                's',NaN);
            %
            %% Solver
            %
            [solverMain,predictorSolver,numJacSolver,defaultSolverOutput] = continuation.solver(obj,obj.stepsizeOptions.rateOfContraction);
            obj.solver = struct('defaultOutput',defaultSolverOutput,...
                                'exitflag',-2,...
                                'jacobian',[],...
                                'main',solverMain,...
                                'numJac',numJacSolver,...
                                'output',defaultSolverOutput,...
                                'predictor',predictorSolver,...
                                'temp',[]);
            %
        end

        function updateOpt(obj)
            %% set direction if l0~=lStart || lTarget~=lEnd
            %
            if ~obj.optIsSet.direction
                if obj.optIsSet.l0
                    if obj.optIsSet.lTarget
                        obj.opt.direction = sign(obj.opt.lTarget-obj.opt.l0)*[zeros(size(obj.info.var0));1];
                    else
                        obj.opt.direction = sign(obj.info.lEnd-obj.opt.l0)*[zeros(size(obj.info.var0));1];
                    end
                else
                    if obj.optIsSet.lTarget
                        obj.opt.direction = sign(obj.opt.lTarget-obj.info.lStart)*[zeros(size(obj.info.var0));1];
                    end
                end
            end
            %
            %% set dsTol dependent on corrector method:
            %
            if ~obj.optIsSet.dsTol
                if obj.opt.corrector.sphere
                    obj.opt.dsTol = [0.99,1.01];
                elseif obj.opt.corrector.orthogonal
                    obj.opt.dsTol = [0.99,5];
                elseif obj.opt.corrector.ellipsoid
                    obj.opt.dsTol = [0.24,1.01];
                elseif obj.opt.corrector.ellipsoid2
                    obj.opt.dsTol = [0.000001,1.01];
                elseif obj.opt.corrector.unique
                    obj.opt.dsTol = [0.99,5];
                elseif obj.opt.corrector.paraboloid
                    obj.opt.dsTol = [0.09,1.01];
                else
                    obj.opt.dsTol = [0.5,1.5];
                end
            end
            %
        end
    end
end