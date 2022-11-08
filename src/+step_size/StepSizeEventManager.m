%% path continuation - step_size.StepSizeEventManager
%
%   Institute of Dynamics and Vibration Research
%   Leibniz University Hannover
%   08.11.2022 - Tido Kubatschek
%
classdef StepSizeEventManager < handle
    %EventManager Manages occuring events in step size control

    properties (Access = private)
        all_events;
        any_active = false;
        active_events = [];
        Initial;
        Current;
        changed = false;
        all_params = [];
    end

    methods
        function ds_new = check_events_and_adapt_stepsize(obj,ds_current,value_struct)
            len = numel(obj.all_events);
            is_active = false(len,1);
            %
            % check all events whether anys condition is met
            %
            obj.active_events = [];
            for k = 1:len
                if obj.all_events{k}.check_event(value_struct)
                    %
                    % additionally check counter
                    %
                    if obj.all_events{k}.check_counter()
                        is_active(k) = true;
                        obj.active_events = [obj.active_events, k];
                    end
                else
                    %
                    % if event was active in the last run but isnt now
                    %
                    if obj.all_events{k}.get('last_active')
                        obj.all_events{k}.set('last_active',false);
                    end
                end
            end
            
            % if there are active events adapt stepsize
            if any(is_active)
                %
                % get new stepsizes
                ds_new = obj.set_stepsize(is_active,ds_current);
                %
                % there are active events
                %
                obj.any_active = true;
                %
            else
                %
                % new stepsize is old stepsize
                %
                ds_new.use = ds_current;
                %
                % set ds_max and ds_min to initial values
                %
                ds_new.ds_min = obj.Initial.ds_min;
                ds_new.ds_max = obj.Initial.ds_max;
                %
                % currently no event is active
                obj.any_active = false;
                %
            end
        end

        function add_event(obj,event_condition, needed_params, ds_min, ds_max, var_input)
            if nargin == 5 || isempty(var_input)
                var_input = [];
            end
            if isempty(obj.all_events)
                obj.all_events = {step_size.StepSizeSingleEvent(event_condition, needed_params, ds_min, ds_max, var_input)};
            else
                len = numel(obj.all_events);
                tmp = cell(len+1,1);
                for k = 1:len
                    tmp{k} = obj.all_events{k};
                end
                tmp{end} = step_size.StepSizeSingleEvent(event_condition, needed_params, ds_min, ds_max, var_input);
                obj.all_events = tmp;
            end
        end

        function output = getActiveEvents(obj)
            output = obj.active_events;
        end
        
        function output = getChanged(obj)
            output = obj.changed;
        end

        function output = getNeededVariables(obj)
            if isempty(obj.all_params)
                list = {};
                for k = 1:numel(obj.all_events)
                    params = obj.all_events{k}.get('needed_parameters');
                    list = [list, params];
                end
                obj.all_params = unique(list);
            end
            output = obj.all_params;
        end

        function obj = set_initial(obj,Initial)
            obj.Initial = Initial;
        end

        function obj = set_current(obj,ds_max,ds_min)
            obj.Current.ds_max = ds_max;
            obj.Current.ds_min = ds_min;
        end
    end
    methods (Access = private)
        function new_stepsizes = set_stepsize(obj,is_active,ds_current)
            ds = inf;
            new_stepsizes.ds_max = [];
            new_stepsizes.ds_min = [];
            obj.changed = false;
            % iterate through events
            for k = 1:numel(obj.all_events)
                % if any is active and its stepsize is lower than the
                % previous saved value and it wasn't active in the latest step
                % save its value
                if is_active(k) && ~obj.all_events{k}.get('last_active') && obj.all_events{k}.get('ds_max') < obj.Current.ds_max
                    obj.changed = true;
                    ds = obj.all_events{k}.get('ds_max');
                    obj.all_events{k}.set('last_active',true);

                    % set min and max stepsize
                    new_stepsizes.ds_max = obj.all_events{k}.get('ds_max');
                    new_stepsizes.ds_min = obj.all_events{k}.get('ds_min');
                end
            end
            %
            % check whether found stepsize is bigger than currently used one
            %
            if ds >= ds_current
                new_stepsizes.ds = ds_current;
            else
                new_stepsizes.ds = ds;
            end
        end
    end
end