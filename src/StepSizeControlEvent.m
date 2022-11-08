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
        
        function obj = addEvent(obj,event_condition,needed_parameters,ds_min,ds_max,counter_max)
            arguments
                obj
                event_condition function_handle
                needed_parameters (1,:) cell
                ds_min (1,1) double {mustBePositive}
                ds_max (1,1) double {mustBePositive,mustBeGreaterThan(ds_max,ds_min)}
                counter_max (1,1) double {mustBeIntegerorInf(counter_max),mustBePositive} = [];
            end
            
            len = length(obj.stepsizeevents);
            idx = len + 1;
            obj.stepsizeevents{idx}.condition = event_condition;
            obj.stepsizeevents{idx}.needed_parameters = needed_parameters;
            obj.stepsizeevents{idx}.ds_min = ds_min;
            obj.stepsizeevents{idx}.ds_max = ds_max;
            
            if nargin == 5
                obj.stepsizeevents{idx}.counter_max = counter_max;
            end
        end

        function event_object = getEvents(obj)
            if ~isempty(obj.stepsizeevents)
                event_object = obj.stepsizeevents;
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