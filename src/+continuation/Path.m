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
        JAll
        lAll
        nL
        nVar
        pathInfoValue
        plus = false
        saveAllJacobian
        signDetJAll
        varAll
    end

    properties (Access = public)
        outputFormat
        speedOfContinuation
        xPredictor
    end

    properties (Access = private)
        plusStruct = struct('bifTestValue',[],'J',[],'pathInfoValue',[],...
            'signDetJ',[],'var',[],'l',[],'xPredictor',[]);
    end

    properties (Dependent)
        nAll
        sAll
        xAll
    end

    %% public methods
    methods
        %% constructor
        function obj = Path(nVar,nL,NameValueArgs)
            arguments
                nVar (1,1) double {mustBeInteger,mustBeGreaterThanOrEqual(nVar,1)}
                nL (1,1) double {mustBeInteger,mustBeGreaterThanOrEqual(nL,1)}
                NameValueArgs.saveAllJacobian (1,1) logical = false
            end
            obj.nVar = nVar;
            obj.nL = nL;
            obj.saveAllJacobian = NameValueArgs.saveAllJacobian;
            obj.bifTestValue = zeros(1,0);
            if obj.saveAllJacobian
                obj.JAll = {};
            else
                obj.JAll = cell(1,3); % initial, last, previous
            end
            obj.lAll = zeros(nL,0);
            obj.pathInfoValue = zeros(1,0);
            obj.signDetJAll = zeros(1,0);
            obj.speedOfContinuation = zeros(1,0);
            obj.varAll = zeros(nVar,0);
            obj.xPredictor = zeros(nVar+nL,1);
            obj.resetOutput();
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
                NameValueArgs.pathInfoValue (1,1) double
            end
            validateattributes(var,{'numeric'},{'size',[obj.nVar,1]});
            validateattributes(l,{'numeric'},{'size',[obj.nL,1]});
            %% add/insert
            obj.varAll = [obj.varAll(:,1:(pos-1)),var,obj.varAll(:,pos:end)];
            obj.lAll = [obj.lAll(:,1:(pos-1)),l,obj.lAll(:,pos:end)];
            if ~obj.saveAllJacobian
                obj.JAll{3} = obj.JAll{2};
                obj.JAll{2} = J;
            end
            obj.signDetJAll = [obj.signDetJAll(:,1:(pos-1)),sign(det(J(1:obj.nVar,1:obj.nVar))),obj.signDetJAll(:,pos:end)];
            if isfield(NameValueArgs,'bifTestValue')
                obj.bifTestValue = [obj.bifTestValue(:,1:(pos-1)),NameValueArgs.bifTestValue,obj.bifTestValue(:,pos:end)];
            end
            if isfield(NameValueArgs,'pathInfoValue')
                obj.pathInfoValue = [obj.pathInfoValue(:,1:(pos-1)),NameValueArgs.pathInfoValue,obj.pathInfoValue(:,pos:end)];
            end
            % ... TODO!
%             obj.bifTestValue = zeros(1,0);
%             obj.pathInfoValue = zeros(1,0);
%             obj.speedOfContinuation = zeros(1,0);
        end

        function addPointAtEnd(obj,var,l,J,NameValueArgs)
            %% arguments
            arguments
                obj (1,1) continuation.Path
                var (:,1) double
                l (1,1) double
                J (:,:) double
                NameValueArgs.bifTestValue (1,1) double
                NameValueArgs.pathInfoValue (1,1) double
            end
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
            if obj.saveAllJacobian
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
            else
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
            % todo...
            error('continuation.Path.remove() is not implemented!');
            obj.bifTestValue(:,idxRmv) = [];
            if obj.saveAllJacobian
                obj.JAll(idxRmv) = [];
            else
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
            end
            obj.lAll(:,idxRmv) = [];
            obj.pathInfoValue(:,idxRmv) = [];
            obj.signDetJAll(:,idxRmv) = [];
            obj.varAll(:,idxRmv) = [];
            obj.xPredictor = [];
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
            if obj.saveAllJacobian
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
                if obj.saveAllJacobian
                    obj.JAll{end+1} = obj.plusStruct.J;
                else
                    obj.JAll{3} = obj.JAll{2};
                    obj.JAll{2} = obj.plusStruct.J;
                end
                obj.lAll = [obj.lAll,obj.plusStruct.l];
                obj.pathInfoValue = [obj.pathInfoValue,obj.plusStruct.pathInfoValue];
                obj.signDetJAll = [obj.signDetJAll,obj.plusStruct.signDetJ];
                obj.varAll = [obj.varAll,obj.plusStruct.var];
                obj.xPredictor = obj.plusStruct.xPredictor;
                %% clear plus struct
                obj.clearPlusStruct();
                %% toggle plus
                obj.plus = false;
            else
                if obj.nAll<2
                    error("Path must have at least two entries to use 'plus'.");
                end
                %% fill plus
                obj.plusStruct.bifTestValue = obj.bifTestValue(:,end);
                if obj.saveAllJacobian
                    obj.plusStruct.J = obj.JAll{end};
                else
                    obj.plusStruct.J = obj.JAll{2};
                end
                obj.plusStruct.l = obj.lAll(:,end);
                obj.plusStruct.pathInfoValue = obj.pathInfoValue(:,end);
                obj.plusStruct.signDetJ = obj.signDetJAll(:,end);
                obj.plusStruct.var = obj.varAll(:,end);
                obj.plusStruct.xPredictor = obj.xPredictor;
                %% cut all
                obj.bifTestValue = obj.bifTestValue(:,1:(end-1));
                if obj.saveAllJacobian
                    obj.JAll = obj.JAll(1:(end-1));
                else
                    obj.JAll{2} = obj.JAll{3};
                    obj.JAll{3} = [];
                end
                obj.lAll = obj.lAll(:,1:(end-1));
                obj.pathInfoValue = obj.pathInfoValue(:,1:(end-1));
                obj.signDetJAll = obj.signDetJAll(:,1:(end-1));
                obj.varAll = obj.varAll(:,1:(end-1));
                obj.xPredictor = [];
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
        function clearPlusStruct(obj)
            obj.plusStruct.bifTestValue = [];
            obj.plusStruct.J = [];
            obj.plusStruct.l = [];
            obj.plusStruct.pathInfoValue = [];
            obj.plusStruct.signDetJ = [];
            obj.plusStruct.var = [];
            obj.plusStruct.xPredictor = [];
        end
    end
end