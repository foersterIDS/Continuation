%% path continuation - StepSizeControlEvent
%
%   Institute of Dynamics and Vibration Research
%   Leibniz University Hannover
%   08.11.2022 - Tido Kubatschek
%
classdef StepSizeControlEvent < handle
    %StepSizeControlEvent Creates event object for event based step size
    %control

    properties (Access = private)
        stepsizeevents = [];
    end

    methods
        function obj = StepSizeControlEvent()
            obj.stepsizeevents = cell(0);
        end
        
        function obj = addEvent(obj,eventCondition,neededParameters,dsMin,dsMax,NameValueArgs)
            arguments
                obj
                eventCondition function_handle
                neededParameters (1,:) cell
                dsMin (1,1) double {mustBePositive}
                dsMax (1,1) double {mustBePositive,mustBeGreaterThan(dsMax,dsMin)}
                NameValueArgs.counterMax (1,1) double {mustBeIntegerorInf(NameValueArgs.counterMax),mustBePositive};
                NameValueArgs.dsAfter (1,1) double {mustBePositive}
            end
            
            len = length(obj.stepsizeevents);
            idx = len + 1;
            obj.stepsizeevents{idx}.condition = eventCondition;
            obj.stepsizeevents{idx}.neededParameters = neededParameters;
            obj.stepsizeevents{idx}.dsMin = dsMin;
            obj.stepsizeevents{idx}.dsMax = dsMax;
            
            if isfield(NameValueArgs,'counterMax')
                obj.stepsizeevents{idx}.counterMax = NameValueArgs.counterMax;
            end
            if isfield(NameValueArgs,'dsAfter')
                obj.stepsizeevents{idx}.dsAfter = NameValueArgs.dsAfter;
            end
        end

        function eventObject = getEvents(obj)
            if ~isempty(obj.stepsizeevents)
                eventObject = obj.stepsizeevents;
            else
                error('Please first create StepSizeControlEvent object and add events to it.');
            end
        end
    end
end

function passed = mustBeIntegerorInf(val)
    if isinf(val) || isinteger(val)
        passed = true;
    else
        passed = false;
    end
end