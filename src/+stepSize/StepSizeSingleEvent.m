%% path continuation - stepSize.StepSizeSingleEvent
%
%   Institute of Dynamics and Vibration Research
%   Leibniz University Hannover
%   08.11.2022 - Tido Kubatschek
%
classdef StepSizeSingleEvent < handle
    % SingleEvent Class for a single event

    properties (Access = private)
        eventCondition
        neededParameters
        dsMin {mustBePositive}
        dsMax {mustBePositive}
        counterMax
        dsAfter

        counter = 0
        lastActive = false
    end

    methods
        function obj = StepSizeSingleEvent(eventCondition, neededParameters, dsMin, dsMax, NameValueArgs)
            arguments
                eventCondition
                neededParameters
                dsMin
                dsMax
                NameValueArgs
            end
            if isa(eventCondition, 'function_handle')
                obj.eventCondition = eventCondition;
            else
                error('eventCondtion must be a function handle.')
            end
            
            for k = 1:length(neededParameters)
                if ~ischar(neededParameters{k})
                    error('needeParameters must be a struct of chars.')
                end
            end
            obj.neededParameters = neededParameters;
            
            obj.dsMin = dsMin;
            obj.dsMax = dsMax;
            
            if isfield(NameValueArgs,'counterMax')
                if (round(NameValueArgs.counterMax) == NameValueArgs.counterMax && NameValueArgs.counterMax >= 0) || isinf(NameValueArgs.counterMax)
                    obj.counterMax = NameValueArgs.counterMax;
                else
                    error('Wrong input in SingleEvent!');
                end
            else
                obj.counterMax = inf;
            end

            if isfield(NameValueArgs,'dsAfter')
                obj.dsAfter = NameValueArgs.dsAfter;
            else
                obj.dsAfter = [];
            end
            
        end

        function output = checkEvent(obj,valueStruct)
            numOfVariables = numel(obj.neededParameters);
            evaluationString = 'obj.eventCondition(';
            for k = 1:numOfVariables
                eval([obj.neededParameters{k}, '= valueStruct.(obj.neededParameters{k})',';']);
                evaluationString = [evaluationString, obj.neededParameters{k}, ','];
            end
            evaluationString(end) = ')';
            output = eval(evaluationString);
        end

        function output = checkCounter(obj)
            output = (obj.counter < obj.counterMax);
        end

        function increaseCounter(obj)
            obj.counter = obj.counter + 1;
        end

        function value = get(obj,propertyName)
            value = obj.(propertyName);
        end

        function obj = set(obj,propertyName,value)
            obj.(propertyName) = value;
        end
    end
end