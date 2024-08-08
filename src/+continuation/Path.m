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
        signDetJAll
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
            switch obj.Opt.jacobianOut
                case 'basic'
                    obj.JAll = cell(1,3); % initial, last, previous
                case 'full'
                    obj.JAll = {};
                otherwise
                    error('jacobianOut must be full or basic!');
            end
            obj.lAll = zeros(nL,0);
            obj.signDetJAll = zeros(1,0);
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
        function addPoint(obj,var,l,J,pos,NameValueArgs)
            %% arguments
            arguments
                obj (1,1) continuation.Path
                var (:,1) double
                l (1,1) double
                J (:,:) double
                pos (1,1) double {mustBeInteger,mustBeGreaterThan(pos,0)} = obj.nAll+1
                NameValueArgs.bifTestValue (1,1) double
                NameValueArgs.iterations (1,1) double {mustBeInteger,mustBeGreaterThan(NameValueArgs.iterations,0)}
                NameValueArgs.pathInfoValue (1,1) double
                NameValueArgs.predictor (:,1) double
                NameValueArgs.rateOfContraction (1,1) double
                NameValueArgs.speedOfContinuation (1,1) double
            end
            validateattributes(var,{'numeric'},{'size',[obj.nVar,1]});
            validateattributes(l,{'numeric'},{'size',[obj.nL,1]});
            obj.checkInputOptions(NameValueArgs);
            %% add/insert
            obj.varAll = [obj.varAll(:,1:(pos-1)),var,obj.varAll(:,pos:end)];
            obj.lAll = [obj.lAll(:,1:(pos-1)),l,obj.lAll(:,pos:end)];
            switch obj.Opt.jacobianOut
                case 'basic'
                    obj.JAll{3} = obj.JAll{2};
                    obj.JAll{2} = J;
                case 'full'
                    obj.JAll = [obj.JAll(:,1:(pos-1)),J,obj.JAll(:,pos:end)];
                otherwise
                    error('jacobianOut must be full or basic!');
            end
            obj.signDetJAll = [obj.signDetJAll(:,1:(pos-1)),sign(det(J(1:obj.nVar,1:obj.nVar))),obj.signDetJAll(:,pos:end)];
            if isfield(NameValueArgs,'bifTestValue')
                obj.bifTestValue = [obj.bifTestValue(:,1:(pos-1)),NameValueArgs.bifTestValue,obj.bifTestValue(:,pos:end)];
            end
            if isfield(NameValueArgs,'iterations')
                obj.iterations = [obj.iterations(:,1:(pos-1)),NameValueArgs.iterations,obj.iterations(:,pos:end)];
            end
            if isfield(NameValueArgs,'pathInfoValue')
                obj.pathInfoValue = [obj.pathInfoValue(:,1:(pos-1)),NameValueArgs.pathInfoValue,obj.pathInfoValue(:,pos:end)];
            end
            if isfield(NameValueArgs,'predictor')
                obj.xPredictorAll = [obj.xPredictorAll(:,1:(pos-1)),NameValueArgs.predictor,obj.xPredictorAll(:,pos:end)];
            end
            if isfield(NameValueArgs,'rateOfContraction')
                obj.rateOfContraction = [obj.rateOfContraction(:,1:(pos-1)),NameValueArgs.rateOfContraction,obj.rateOfContraction(:,pos:end)];
            end
            if isfield(NameValueArgs,'speedOfContinuation')
                obj.speedOfContinuation = [obj.speedOfContinuation(:,1:(pos-1)),NameValueArgs.speedOfContinuation,obj.speedOfContinuation(:,pos:end)];
            end
        end

        function addPointAtEnd(obj,var,l,J,NameValueArgs)
            %% arguments
            arguments
                obj (1,1) continuation.Path
                var (:,1) double
                l (1,1) double
                J (:,:) double
                NameValueArgs.bifTestValue (1,1) double
                NameValueArgs.iterations (1,1) double {mustBeInteger,mustBeGreaterThan(NameValueArgs.iterations,0)}
                NameValueArgs.pathInfoValue (1,1) double
                NameValueArgs.predictor (:,1) double
                NameValueArgs.rateOfContraction (1,1) double
                NameValueArgs.speedOfContinuation (1,1) double
            end
            obj.checkInputOptions(NameValueArgs);
            %% pass to addPoint(...)
            pos = obj.nAll+1;
            nva = aux.struct2cellPreserveFieldnames(NameValueArgs);
            obj.addPoint(var,l,J,pos,nva{:});
        end

        function J = getJacobianByName(obj,jacName)
            %% arguments
            arguments
                obj (1,1) continuation.Path
                jacName (1,:) char {mustBeMember(jacName,{'last','previous','initial','plus'})}
            end
            %% get jacobian
            switch obj.Opt.jacobianOut
                case 'basic'
                    switch jacName
                        case 'initial'
                            J = obj.JAll{1};
                        case 'last'
                            J = obj.JAll{2};
                        case 'plus'
                            J = obj.plusStruct.J;
                        case 'previous'
                            J = obj.JAll{3};
                    end
                case 'full'
                    switch jacName
                        case 'initial'
                            J = obj.JAll{1};
                        case 'last'
                            J = obj.JAll{end};
                        case 'plus'
                            J = obj.plusStruct.J;
                        case 'previous'
                            J = obj.JAll{end-1};
                    end
                otherwise
                    error('jacobianOut must be full or basic!');
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
            switch obj.Opt.jacobianOut
                case 'basic'
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
                case 'full'
                    obj.JAll(idxRmv) = [];
                otherwise
                    error('jacobianOut must be full or basic!');
            end
            obj.lAll(:,idxRmv) = [];
            obj.signDetJAll(:,idxRmv) = [];
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
            switch obj.Opt.jacobianOut
                case 'basic'
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
                case 'full'
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
                otherwise
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
                switch obj.Opt.jacobianOut
                    case 'basic'
                        obj.JAll{3} = obj.JAll{2};
                        obj.JAll{2} = obj.plusStruct.J;
                    case 'full'
                        obj.JAll{end+1} = obj.plusStruct.J;
                    otherwise
                        error('jacobianOut must be full or basic!');
                end
                obj.lAll = [obj.lAll,obj.plusStruct.l];
                obj.pathInfoValue = [obj.pathInfoValue,obj.plusStruct.pathInfoValue];
                obj.signDetJAll = [obj.signDetJAll,obj.plusStruct.signDetJ];
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
                switch obj.Opt.jacobianOut
                    case 'basic'
                        obj.plusStruct.J = obj.JAll{2};
                    case 'full'
                        obj.plusStruct.J = obj.JAll{end};
                    otherwise
                        error('jacobianOut must be full or basic!');
                end
                obj.plusStruct.l = obj.lAll(:,end);
                obj.plusStruct.signDetJ = obj.signDetJAll(:,end);
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
                switch obj.Opt.jacobianOut
                    case 'basic'
                        obj.JAll{2} = obj.JAll{3};
                        obj.JAll{3} = [];
                    case 'full'
                        obj.JAll = obj.JAll(1:(end-1));
                    otherwise
                        error('jacobianOut must be full or basic!');
                end
                obj.lAll = obj.lAll(:,1:(end-1));
                obj.signDetJAll = obj.signDetJAll(:,1:(end-1));
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
            if obj.StepsizeOptions.iterations && ~isfield(nvaStruct,'iterations')
                error('A iterations must be passed if a iterations is set.')
            end
            if obj.StepsizeOptions.predictor && ~isfield(nvaStruct,'predictor')
                error('A predictor must be passed if a predictor is set.')
            end
            if obj.StepsizeOptions.rateOfContraction && ~isfield(nvaStruct,'rateOfContraction')
                error('A rateOfContraction must be passed if a rateOfContraction is set.')
            end
            if obj.StepsizeOptions.speedOfContinuation && ~isfield(nvaStruct,'speedOfContinuation')
                error('A speedOfContinuation must be passed if a speedOfContinuation is set.')
            end
        end

        function clearPlusStruct(obj)
            obj.plusStruct.bifTestValue = [];
            obj.plusStruct.iterations = [];
            obj.plusStruct.J = [];
            obj.plusStruct.l = [];
            obj.plusStruct.pathInfoValue = [];
            obj.plusStruct.rateOfContraction = [];
            obj.plusStruct.signDetJ = [];
            obj.plusStruct.speedOfContinuation = [];
            obj.plusStruct.var = [];
            obj.plusStruct.xPredictor = [];
        end
    end
end