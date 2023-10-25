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
        nv
        xAll
    end
    %% methods
    methods
        %% constructor
        function obj = Path(NameValueArgs)
            arguments
                NameValueArgs.saveAllJacobian (1,1) logical = false
            end
            obj.saveAllJacobian = NameValueArgs.saveAllJacobian;
            obj.biftestValue = [];
            if obj.saveAllJacobian
                obj.JAll = {};
            else
                obj.JAll = cell(1,3); % initial, last, previous
            end
            obj.lAll = [];
            obj.pathInfoValue = [];
            obj.sAll = [];
            obj.signDetJAll = [];
            obj.speedOfContinuation = [];
            obj.varAll = [];
            obj.xPredictor = [];
        end
        %% getter
        function nAll = get.nAll(obj)
            nAll = numel(obj.lAll);
        end

        function nv = get.nv(obj)
            nv = numel(obj.varAll(:,1));
        end

        function xAll = get.xAll(obj)
            xAll = [obj.varAll;obj.lAll];
        end
        %% setter
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
    end
end