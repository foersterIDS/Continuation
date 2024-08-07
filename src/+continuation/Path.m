%% path continuation - continuation.Path (class)
%
%   Institute of Dynamics and Vibration Research
%   Leibniz University Hannover
%   28.07.2023 - Alwin Förster
%
classdef Path < handle
    %% properties
    properties (GetAccess = public, SetAccess = private)
        biftestValue
        JAll
        lAll
        nL
        nVar
        pathInfoValue
        plus = false
        sAll
        saveAllJacobian
        signDetJAll
        varAll
    end

    properties (Access = public)
        outputFormat
        speedOfContinuation
        xPredictor
    end

    properties (Dependent)
        nAll
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
            obj.biftestValue = zeros(1,0);
            if obj.saveAllJacobian
                obj.JAll = {};
            else
                obj.JAll = cell(1,3); % initial, last, previous
            end
            obj.lAll = zeros(nL,0);
            obj.pathInfoValue = zeros(1,0);
            obj.sAll = zeros(1,0);
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
                xP (1,:) double {mustBeVector}
            end
            xP = xP(:);
            validateattributes(xP,{'double'},{'size',[obj.nVar+1,1]});
            obj.xPredictor = xP;
        end

        %% general methods
        function togglePlus(obj)
            %% arguments
            arguments
                obj (1,1) continuation.Path
            end
            %% toggle plus
            if plus
                % todo...
                % xPlus löschen und wieder hinten an Pfad anfügen
            else
                % todo...
                % xPlus befüllen und von Pfad entfernen
            end
        end

        function addPoint(obj,var,l,J,pos,NameValueArgs)
            %% arguments
            arguments
                obj (1,1) continuation.Path
                var (:,1) double
                l (1,1) double
                J (:,:) double
                pos (1,1) double {mustBeInteger,mustBeGreaterThan(pos,0)} = obj.nAll+1
                NameValueArgs.biftestValue (1,1) double
            end
            validateattributes(var,{'numeric'},{'size',[obj.nVar,1]});
            validateattributes(l,{'numeric'},{'size',[obj.nL,1]});
            %% add/insert
            obj.varAll = [obj.varAll(:,1:(pos-1)),var,obj.varAll(:,pos:end)];
            obj.lAll = [obj.lAll(:,1:(pos-1)),l,obj.lAll(:,pos:end)];
            obj.sAll = [0,cumsum(sqrt(sum(diff(obj.xAll,1,2).^2,1)))];
            if ~obj.saveAllJacobian
                obj.JAll{3} = obj.JAll{2};
                obj.JAll{2} = J;
            end
            obj.signDetJAll = [obj.signDetJAll(:,1:(pos-1)),sign(det(J(1:obj.nVar,1:obj.nVar))),obj.signDetJAll(:,pos:end)];
            if isfield(NameValueArgs,'biftestValue')
                obj.biftestValue = [obj.biftestValue(:,1:(pos-1)),NameValueArgs.biftestValue,obj.biftestValue(:,pos:end)];
            end
            % ... TODO!
%             obj.biftestValue = zeros(1,0);
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
                NameValueArgs.biftestValue (1,1) double
            end
            %% pass to addPoint(...)
            pos = obj.nAll+1;
            nva = aux.struct2cellPreserveFieldnames(NameValueArgs);
            obj.addPoint(var,l,J,pos,nva{:});
        end

        function li = lIndex(obj,index)
            arguments
                obj (1,1) continuation.Path
                index (1,:) double {mustBeInteger,mustBeGreaterThanOrEqual(index,1)}
            end
            validateattributes(index,{'double'},'<=',obj.nAll);
            li = obj.lAll(:,index);
        end

        function resetOutput(obj)
            obj.outputFormat = 'singleValueForL';
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
        % ...
    end
end