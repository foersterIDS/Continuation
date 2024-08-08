%% path continuation - continuation.Path (class)
%
%   Institute of Dynamics and Vibration Research
%   Leibniz University Hannover
%   28.07.2023 - Alwin Förster
%
classdef Path < handle
    %% properties
    properties (GetAccess = public, SetAccess = private)
        bifTestValue
        iterations
        JAll
        lAll
        nL
        nVar
        pathInfoValue
        plus = false
        rateOfContraction
        signDetJRedAll
        speedOfContinuation
        varAll
        xPredictorAll
    end

    properties (Access = public)
        outputFormat
        xPredictor
    end

    properties (Access = private)
        Opt
        OptIsSet
        plusStruct = struct('bifTestValue',[],'iterations',[],'J',[],...
                            'pathInfoValue',[],'rateOfContraction',[],...
                            'signDetJ',[],'speedOfContinuation',[],...
                            'var',[],'l',[],'xPredictor',[]);
        StepsizeOptions
    end

    properties (Dependent)
        nAll
        sAll
        xAll
    end

    %% public methods
    methods
        %% constructor
        function obj = Path(nVar,nL,Opt,OptIsSet,StepsizeOptions)
            %% arguments
            arguments
                nVar (1,1) double {mustBeInteger,mustBeGreaterThanOrEqual(nVar,1)}
                nL (1,1) double {mustBeInteger,mustBeGreaterThanOrEqual(nL,1)}
                Opt (1,1) struct
                OptIsSet (1,1) struct
                StepsizeOptions (1,1) struct
            end
            %% initialize
            obj.nVar = nVar;
            obj.nL = nL;
            obj.Opt = Opt;
            obj.OptIsSet = OptIsSet;
            obj.StepsizeOptions = StepsizeOptions;
            if obj.Opt.jacobianOut.basic
                obj.JAll = cell(1,3); % initial, last, previous
            elseif obj.Opt.jacobianOut.full
                obj.JAll = {};
            else
                error('jacobianOut must be full or basic!');
            end
            obj.lAll = zeros(nL,0);
            obj.signDetJRedAll = zeros(1,0);
            obj.varAll = zeros(nVar,0);
            obj.xPredictor = zeros(nVar+nL,1);
            obj.resetOutput();
            %% additional functions
            if obj.OptIsSet.pathInfoFunction
                obj.pathInfoValue = zeros(1,0);
            end
            if obj.OptIsSet.bifAdditionalTestfunction
                obj.bifTestValue = zeros(1,0);
            end
            %% step size control values
            if obj.StepsizeOptions.iterations
                obj.iterations = zeros(1,0);
            end
            if obj.StepsizeOptions.predictor
                obj.xPredictorAll = zeros(nVar+nL,0);
            end
            if obj.StepsizeOptions.rateOfContraction
                obj.rateOfContraction = zeros(1,0);
            end
            if obj.StepsizeOptions.speedOfContinuation
                obj.speedOfContinuation = zeros(1,0);
            end
        end

        %% getter
        function nAll = get.nAll(obj)
            nAll = numel(obj.lAll(1,:));
        end

        function sAll = get.sAll(obj)
            sAll = [0,cumsum(sqrt(sum(diff(obj.xAll,1,2).^2,1)))];
        end

        function xAll = get.xAll(obj)
            xAll = [obj.varAll;obj.lAll];
        end

        %% setter
        function set.outputFormat(obj,oF)
            arguments
                obj (1,1) continuation.Path
                oF (1,:) char {mustBeMember(oF,{'singleValueForL','full'})}
            end
            obj.outputFormat = oF;
        end

        function set.xPredictor(obj,xP)
            arguments
                obj (1,1) continuation.Path
                xP (:,1) double {mustBeVector}
            end
            xP = xP(:);
            % validateattributes(xP,{'double'},{'size',[obj.nVar+1,1]});
            obj.xPredictor = xP;
        end

        %% general methods
        function addPoint(obj,var,l,J,Solver,pos,NameValueArgs)
            %% arguments
            arguments
                obj (1,1) continuation.Path
                var (:,1) double
                l (1,1) double
                J (:,:) double
                Solver (1,1) struct
                pos (1,1) double {mustBeInteger,mustBeGreaterThan(pos,0)} = obj.nAll+1
                NameValueArgs.bifTestValue (1,1) double
                NameValueArgs.pathInfoValue (1,1) double
                NameValueArgs.predictor (:,1) double
                NameValueArgs.speedOfContinuation (1,1) double
            end
            validateattributes(var,{'numeric'},{'size',[obj.nVar,1]});
            validateattributes(l,{'numeric'},{'size',[obj.nL,1]});
            obj.checkInputOptions(NameValueArgs);
            %% add/insert
            obj.varAll = [obj.varAll(:,1:(pos-1)),var,obj.varAll(:,pos:end)];
            obj.lAll = [obj.lAll(:,1:(pos-1)),l,obj.lAll(:,pos:end)];
            if obj.Opt.jacobianOut.basic
                obj.JAll{3} = obj.JAll{2};
                obj.JAll{2} = J;
            elseif obj.Opt.jacobianOut.full
                obj.JAll = [obj.JAll(:,1:(pos-1)),J,obj.JAll(:,pos:end)];
            else
                error('jacobianOut must be full or basic!');
            end
            obj.signDetJRedAll = [obj.signDetJRedAll(:,1:(pos-1)),sign(det(J(1:obj.nVar,1:obj.nVar))),obj.signDetJRedAll(:,pos:end)];
            if obj.OptIsSet.bifAdditionalTestfunction
                obj.bifTestValue = [obj.bifTestValue(:,1:(pos-1)),NameValueArgs.bifTestValue,obj.bifTestValue(:,pos:end)];
            end
            if obj.StepsizeOptions.iterations
                obj.iterations = [obj.iterations(:,1:(pos-1)),Solver.output.iterations,obj.iterations(:,pos:end)];
            end
            if obj.OptIsSet.pathInfoFunction
                obj.pathInfoValue = [obj.pathInfoValue(:,1:(pos-1)),NameValueArgs.pathInfoValue,obj.pathInfoValue(:,pos:end)];
            end
            if obj.StepsizeOptions.predictor
                obj.xPredictorAll = [obj.xPredictorAll(:,1:(pos-1)),NameValueArgs.predictor,obj.xPredictorAll(:,pos:end)];
            end
            if obj.StepsizeOptions.rateOfContraction
                obj.rateOfContraction = [obj.rateOfContraction(:,1:(pos-1)),Solver.output.rateOfContraction,obj.rateOfContraction(:,pos:end)];
            end
            if obj.StepsizeOptions.speedOfContinuation
                obj.speedOfContinuation = [obj.speedOfContinuation(:,1:(pos-1)),NameValueArgs.speedOfContinuation,obj.speedOfContinuation(:,pos:end)];
            end
        end

        function addPointAtEnd(obj,var,l,J,Solver,NameValueArgs)
            %% arguments
            arguments
                obj (1,1) continuation.Path
                var (:,1) double
                l (1,1) double
                J (:,:) double
                Solver (1,1) struct
                NameValueArgs.bifTestValue (1,1) double
                NameValueArgs.pathInfoValue (1,1) double
                NameValueArgs.predictor (:,1) double
                NameValueArgs.speedOfContinuation (1,1) double
            end
            obj.checkInputOptions(NameValueArgs);
            %% pass to addPoint(...)
            pos = obj.nAll+1;
            nva = aux.struct2cellPreserveFieldnames(NameValueArgs);
            obj.addPoint(var,l,J,Solver,pos,nva{:});
        end

        function J = getJacobianByName(obj,jacName,idx1,idx2)
            %% arguments
            arguments
                obj (1,1) continuation.Path
                jacName (1,:) char {mustBeMember(jacName,{'last','previous','initial','plus'})}
                idx1 (1,:) double {mustBeInteger,mustBeGreaterThan(idx1,0)} = 1:obj.nVar
                idx2 (1,:) double {mustBeInteger,mustBeGreaterThan(idx2,0)} = 1:(obj.nVar+obj.nL)
            end
            %% get jacobian
            if obj.Opt.jacobianOut.basic
                switch jacName
                    case 'initial'
                        J = obj.JAll{1};
                    case 'last'
                        J = obj.JAll{2};
                    case 'plus'
                        J = obj.plusStruct.J;
                    case 'previous'
                        if obj.nAll==1
                            J = obj.JAll{2};
                        else
                            J = obj.JAll{3};
                        end
                end
            elseif obj.Opt.jacobianOut.full
                switch jacName
                    case 'initial'
                        J = obj.JAll{1};
                    case 'last'
                        J = obj.JAll{end};
                    case 'plus'
                        J = obj.plusStruct.J;
                    case 'previous'
                        if obj.nAll==1
                            J = obj.JAll{end};
                        else
                            J = obj.JAll{end-1};
                        end
                end
            else
                error('jacobianOut must be full or basic!');
            end
            if ~isinf(idx1)
                idx1(logical((idx1>size(J,1))+(idx1<1))) = [];
                J = J(idx1,:);
            end
            if ~isinf(idx2)
                idx2(logical((idx2>size(J,2))+(idx2<1))) = [];
                J = J(:,idx2);
            end
        end

        function li = lIndex(obj,index)
            arguments
                obj (1,1) continuation.Path
                index (1,:) double {mustBeInteger,mustBeGreaterThanOrEqual(index,1)}
            end
            validateattributes(index,{'double'},'<=',obj.nAll);
            li = obj.lAll(:,index);
        end

        function remove(obj,idxRmv)
            %% arguments
            arguments
                obj (1,1) continuation.Path
                idxRmv (1,:) double {mustBeInteger,mustBeGreaterThan(idxRmv,0)}
            end
            %% remove
            if obj.Opt.jacobianOut.basic
                if sum(idxRmv==obj.nAll)
                    if sum(idxRmv==obj.nAll-1)
                        obj.JAll{2} = [];
                        obj.JAll{3} = [];
                    else
                        obj.JAll{2} = obj.JAll{3};
                        obj.JAll{3} = [];
                    end
                elseif sum(idxRmv==obj.nAll-1)
                    obj.JAll{3} = [];
                end
                if sum(idxRmv==1)
                    obj.JAll{1} = [];
                end
            elseif obj.Opt.jacobianOut.full
                obj.JAll(idxRmv) = [];
            else
                error('jacobianOut must be full or basic!');
            end
            obj.lAll(:,idxRmv) = [];
            obj.signDetJRedAll(:,idxRmv) = [];
            obj.varAll(:,idxRmv) = [];
            obj.xPredictor = [];
            if obj.OptIsSet.bifAdditionalTestfunction
                obj.bifTestValue(:,idxRmv) = [];
            end
            if obj.StepsizeOptions.iterations
                obj.iterations(:,idxRmv) = [];
            end
            if obj.OptIsSet.pathInfoFunction
                obj.pathInfoValue(:,idxRmv) = [];
            end
            if obj.StepsizeOptions.predictor
                obj.xPredictorAll(:,idxRmv) = [];
            end
            if obj.StepsizeOptions.rateOfContraction
                obj.rateOfContraction(:,idxRmv) = [];
            end
            if obj.StepsizeOptions.speedOfContinuation
                obj.speedOfContinuation(:,idxRmv) = [];
            end
        end

        function resetOutput(obj)
            obj.outputFormat = 'singleValueForL';
        end

        function setJacobianByName(obj,J,jacName)
            %% arguments
            arguments
                obj (1,1) continuation.Path
                J (:,:) double
                jacName (1,:) char {mustBeMember(jacName,{'last','previous','initial','plus'})}
            end
            %% get jacobian
            if obj.Opt.jacobianOut.basic
                switch jacName
                    case 'initial'
                        obj.JAll{1} = J;
                    case 'last'
                        obj.JAll{2} = J;
                    case 'plus'
                        obj.plusStruct.J = J;
                    case 'previous'
                        obj.JAll{3} = J;
                end
            elseif obj.Opt.jacobianOut.full
                switch jacName
                    case 'initial'
                        obj.JAll{1} = J;
                    case 'last'
                        obj.JAll{end} = J;
                    case 'plus'
                        obj.plusStruct.J = J;
                    case 'previous'
                        obj.JAll{end-1} = J;
                end
            else
                error('jacobianOut must be full or basic!');
            end
        end

        function suspend(obj)
            obj.clearPlusStruct();
        end

        function toggleStepback(obj)
            %% arguments
            arguments
                obj (1,1) continuation.Path
            end
            %% toggle plus
            if obj.plus
                %% add plus
                obj.bifTestValue = [obj.bifTestValue,obj.plusStruct.bifTestValue];
                if obj.Opt.jacobianOut.basic
                    obj.JAll{3} = obj.JAll{2};
                    obj.JAll{2} = obj.plusStruct.J;
                elseif obj.Opt.jacobianOut.full
                    obj.JAll{end+1} = obj.plusStruct.J;
                else
                    error('jacobianOut must be full or basic!');
                end
                obj.lAll = [obj.lAll,obj.plusStruct.l];
                obj.pathInfoValue = [obj.pathInfoValue,obj.plusStruct.pathInfoValue];
                obj.signDetJRedAll = [obj.signDetJRedAll,obj.plusStruct.signDetJRed];
                obj.varAll = [obj.varAll,obj.plusStruct.var];
                obj.xPredictor = obj.plusStruct.xPredictor;
                if obj.OptIsSet.bifAdditionalTestfunction
                    obj.bifTestValue = [obj.bifTestValue,obj.plusStruct.bifTestValue];
                end
                if obj.StepsizeOptions.iterations
                    obj.iterations = [obj.iterations,obj.plusStruct.iterations];
                end
                if obj.OptIsSet.pathInfoFunction
                    obj.pathInfoValue = [obj.pathInfoValue,obj.plusStruct.pathInfoValue];
                end
                if obj.StepsizeOptions.predictor
                    obj.xPredictorAll = [obj.xPredictorAll,obj.plusStruct.xPredictor];
                end
                if obj.StepsizeOptions.rateOfContraction
                    obj.rateOfContraction = [obj.rateOfContraction,obj.plusStruct.rateOfContraction];
                end
                if obj.StepsizeOptions.speedOfContinuation
                    obj.speedOfContinuation = [obj.speedOfContinuation,obj.plusStruct.speedOfContinuation];
                end
                %% clear plus struct
                obj.clearPlusStruct();
                %% toggle plus
                obj.plus = false;
            else
                if obj.nAll<2
                    error("Path must have at least two entries to use 'plus'.");
                end
                %% fill plus
                if obj.Opt.jacobianOut.basic
                    obj.plusStruct.J = obj.JAll{2};
                elseif obj.Opt.jacobianOut.full
                    obj.plusStruct.J = obj.JAll{end};
                else
                    error('jacobianOut must be full or basic!');
                end
                obj.plusStruct.l = obj.lAll(:,end);
                obj.plusStruct.signDetJRed = obj.signDetJRedAll(:,end);
                obj.plusStruct.var = obj.varAll(:,end);
                if obj.OptIsSet.bifAdditionalTestfunction
                    obj.plusStruct.bifTestValue = obj.bifTestValue(:,end);
                end
                if obj.StepsizeOptions.iterations
                    obj.plusStruct.iterations = obj.iterations(:,end);
                end
                if obj.OptIsSet.pathInfoFunction
                    obj.plusStruct.pathInfoValue = obj.pathInfoValue(:,end);
                end
                if obj.StepsizeOptions.predictor
                    obj.plusStruct.xPredictor = obj.xPredictorAll(:,end);
                end
                if obj.StepsizeOptions.rateOfContraction
                    obj.plusStruct.rateOfContraction = obj.rateOfContraction(:,end);
                end
                if obj.StepsizeOptions.speedOfContinuation
                    obj.plusStruct.speedOfContinuation = obj.speedOfContinuation(:,end);
                end
                %% cut all
                if obj.Opt.jacobianOut.basic
                    obj.JAll{2} = obj.JAll{3};
                    obj.JAll{3} = [];
                elseif obj.Opt.jacobianOut.full
                    obj.JAll = obj.JAll(1:(end-1));
                else
                    error('jacobianOut must be full or basic!');
                end
                obj.lAll = obj.lAll(:,1:(end-1));
                obj.signDetJRedAll = obj.signDetJRedAll(:,1:(end-1));
                obj.varAll = obj.varAll(:,1:(end-1));
                obj.xPredictor = [];
                if obj.OptIsSet.bifAdditionalTestfunction
                    obj.bifTestValue = obj.bifTestValue(:,1:(end-1));
                end
                if obj.StepsizeOptions.iterations
                    obj.iterations = obj.iterations(:,1:(end-1));
                end
                if obj.OptIsSet.pathInfoFunction
                    obj.pathInfoValue = obj.pathInfoValue(:,1:(end-1));
                end
                if obj.StepsizeOptions.predictor
                    obj.xPredictorAll = obj.xPredictorAll(:,1:(end-1));
                end
                if obj.StepsizeOptions.rateOfContraction
                    obj.rateOfContraction = obj.rateOfContraction(:,1:(end-1));
                end
                if obj.StepsizeOptions.speedOfContinuation
                    obj.speedOfContinuation = obj.speedOfContinuation(:,1:(end-1));
                end
                %% toggle plus
                obj.plus = true;
            end
        end

        function turn(obj)
            %% arguments
            arguments
                obj (1,1) continuation.Path
            end
            %% deactivate stepback
            obj.plus = false;
            obj.clearPlusStruct();
            %% turn path
            idxTurn = obj.nAll:-1:1;
            if obj.Opt.jacobianOut.basic
                JLastTemp = obj.JAll{2};
                obj.JAll{2} = obj.JAll{1};
                obj.JAll{3} = [];
                obj.JAll{1} = JLastTemp;
            elseif obj.Opt.jacobianOut.full
                obj.JAll = obj.JAll(idxTurn);
            else
                error('jacobianOut must be full or basic!');
            end
            obj.lAll = obj.lAll(:,idxTurn);
            obj.signDetJRedAll = obj.signDetJRedAll(:,idxTurn);
            obj.varAll = obj.varAll(:,idxTurn);
            obj.xPredictor = [];
            if obj.OptIsSet.bifAdditionalTestfunction
                obj.bifTestValue = obj.bifTestValue(:,idxTurn);
            end
            if obj.StepsizeOptions.iterations
                obj.iterations = obj.iterations(:,idxTurn);
            end
            if obj.OptIsSet.pathInfoFunction
                obj.pathInfoValue = obj.pathInfoValue(:,idxTurn);
            end
            if obj.StepsizeOptions.predictor
                obj.xPredictorAll = obj.xPredictorAll(:,idxTurn);
            end
            if obj.StepsizeOptions.rateOfContraction
                obj.rateOfContraction = obj.rateOfContraction(:,idxTurn);
            end
            if obj.StepsizeOptions.speedOfContinuation
                obj.speedOfContinuation = obj.speedOfContinuation(:,idxTurn);
            end
        end

        function vari = varIndex(obj,index)
            arguments
                obj (1,1) continuation.Path
                index (1,:) double {mustBeInteger,mustBeGreaterThanOrEqual(index,1)}
            end
            validateattributes(index,{'double'},'<=',obj.nAll);
            vari = obj.varAll(:,index);
        end

        function xi = xIndex(obj,index)
            arguments
                obj (1,1) continuation.Path
                index (1,:) double {mustBeInteger,mustBeGreaterThanOrEqual(index,1)}
            end
            validateattributes(index,{'double'},'<=',obj.nAll);
            xi = obj.xAll(:,index);
        end
    end

    %% private methods
    methods (Access = private)
        function checkInputOptions(obj,nvaStruct)
            if obj.OptIsSet.pathInfoFunction && ~isfield(nvaStruct,'pathInfoValue')
                error('A pathInfoValue must be passed if a pathInfoFunction is set.')
            end
            if obj.OptIsSet.bifAdditionalTestfunction && ~isfield(nvaStruct,'bifTestValue')
                error('A bifTestValue must be passed if a bifAdditionalTestfunction is set.')
            end
            if obj.StepsizeOptions.predictor && ~isfield(nvaStruct,'predictor')
                error('predictor must be passed if predictor is set.')
            end
            if obj.StepsizeOptions.speedOfContinuation && ~isfield(nvaStruct,'speedOfContinuation')
                error('speedOfContinuation must be passed if speedOfContinuation is set.')
            end
        end

        function clearPlusStruct(obj)
            obj.plusStruct.bifTestValue = [];
            obj.plusStruct.iterations = [];
            obj.plusStruct.J = [];
            obj.plusStruct.l = [];
            obj.plusStruct.pathInfoValue = [];
            obj.plusStruct.rateOfContraction = [];
            obj.plusStruct.signDetJRed = [];
            obj.plusStruct.speedOfContinuation = [];
            obj.plusStruct.var = [];
            obj.plusStruct.xPredictor = [];
        end
    end
end