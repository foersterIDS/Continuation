%% path continuation - step_size.StepSizeSingleEvent
%
%   Institute of Dynamics and Vibration Research
%   Leibniz University Hannover
%   08.11.2022 - Tido Kubatschek
%
classdef StepSizeSingleEvent < handle
    % SingleEvent Class for a single event

    properties (Access = private)
        event_condition
        needed_parameters
        ds_min {mustBePositive}
        ds_max {mustBePositive}
        counter_max

        counter = 0
        last_active = false
    end

    methods
        function obj = StepSizeSingleEvent(event_condition, needed_parameters, ds_min, ds_max, variable_input)
            if isa(event_condition, 'function_handle')
                obj.event_condition = event_condition;
            else
                error('event_condtion must be a function handle.')
            end
            
            for k = 1:length(needed_parameters)
                if ~ischar(needed_parameters{k})
                    error('neede_parameters must be a struct of chars.')
                end
            end
            obj.needed_parameters = needed_parameters;
            
            obj.ds_min = ds_min;
            obj.ds_max = ds_max;
            
            if ~isempty(variable_input)
                if strcmp(variable_input{1},'counter')
                    if (round(variable_input{2}) == variable_input{2} && variable_input{2} >= 0) || isinf(variable_input{2})
                        obj.counter_max = variable_input{2};
                    else
                        error('Wrong input in SingleEvent!');
                    end
                else
                    error('Wrong input in SingleEvent!');
                end
            else
                obj.counter_max = inf;
            end
            
        end

        function output = check_event(obj,value_struct)
            num_of_variables = numel(obj.needed_parameters);
            evaluation_string = 'obj.event_condition(';
            for k = 1:num_of_variables
                eval([obj.needed_parameters{k}, '= value_struct.(obj.needed_parameters{k})',';']);
                evaluation_string = [evaluation_string, obj.needed_parameters{k}, ','];
            end
            evaluation_string(end) = ')';
            output = eval(evaluation_string);
        end

        function output = check_counter(obj)
            output = (obj.counter < obj.counter_max);
        end

        function value = get(obj,property_name)
            value = obj.(property_name);
        end

        function obj = set(obj,property_name,value)
            obj.(property_name) = value;
        end
    end
end