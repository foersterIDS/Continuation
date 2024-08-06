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
        outputFormat
        pathInfoValue
        sAll
        saveAllJacobian
        signDetJAll
        speedOfContinuation
        varAll
        xPredictor
    end

    properties (Dependent)
        nAll
        xAll
    end
    %% methods
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
            obj.xPredictor = zeros(nVar+nL,0);
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
                obj (1,1) Path
                oF (1,:) char {mustBeMember(oF,{'singleValueForL','full'})}
            end
            obj.outputFormat = oF;
        end

        %% general methods
        function addPoint(obj,var,l,J,pos,NameValueArgs)
            %% arguments
            arguments
                obj (1,1) Path
                var (:,1) double
                l (1,1) double
                J (:,:) double
                pos (1,1) double {mustBeInteger,mustBeGreaterThan(pos,0)} = obj.nAll+1
                NameValueArgs.biftestValue (1,1) double
            end
            %% init.
            nAll = obj.nAll;
            if pos==nAll+1
                doInsert = false;
            elseif pos>0 && pos<=nAll
                doInsert = true;
            else
                error('pos must be a scalar double between 1 and nAll+1.');
            end
            %% add/insert
            % ... TODO!
        end

        function resetOutput(obj)
            obj.outputFormat = 'singleValueForL';
        end
    end
end