%% path continuation - stepSize.StepSizeEventManager
%
%   Institute of Dynamics and Vibration Research
%   Leibniz University Hannover
%   08.11.2022 - Tido Kubatschek
%
classdef StepSizeEventManager < handle
    %EventManager Manages occuring events in step size control

    properties (Access = private)
        allEvents;
        anyActive = false;
        activeEvents = [];
        Initial;
        Current;
        changed = false;
        allParams = [];
    end

    methods
        function dsNew = checkEventsAndAdaptStepsize(obj,dsCurrent,valueStruct)
            len = numel(obj.allEvents);
            isActive = false(len,1);
            %
            % check all events whether anys condition is met
            %
            obj.activeEvents = [];
            kLastActive = zeros(len,1);
            for k = 1:len
                if obj.allEvents{k}.checkEvent(valueStruct)
                    %
                    % additionally check counter
                    %
                    if obj.allEvents{k}.checkCounter()
                        isActive(k) = true;
                        obj.activeEvents = [obj.activeEvents, k];
                    end
                else
                    %
                    % if event was active in the last run but isnt now
                    %
                    if obj.allEvents{k}.get('lastActive')
                        kLastActive(k) = 1;
                        obj.allEvents{k}.set('lastActive',false);
                    end
                end
            end
            
            % if there are active events adapt stepsize
            if any(isActive)
                %
                % get new stepsizes
                dsNew = obj.setStepsize(isActive,dsCurrent);
                %
                % there are active events
                %
                obj.anyActive = true;
                %
            else
                %
                % new stepsize is old stepsize
                %
                obj.changed = false;
                if any(kLastActive)
                    for kk = 1:len
                        if kLastActive(kk) && ~isempty(obj.allEvents{kk}.get('dsAfter'))
                            dsNew.ds = obj.allEvents{kk}.get('dsAfter');
                            obj.changed = true;
                        else
                            dsNew.ds = dsCurrent;
                        end
                    end
                else
                    dsNew.ds = dsCurrent;
                end
                %
                % set dsMax and dsMin to initial values
                %
                dsNew.dsMin = obj.Initial.dsMin;
                dsNew.dsMax = obj.Initial.dsMax;
                %
                % currently no event is active
                obj.anyActive = false;
                %
            end
        end

        function addEvent(obj,eventCondition, neededParams, dsMin, dsMax, NameValueArgs)
            arguments
                obj
                eventCondition
                neededParams
                dsMin
                dsMax
                NameValueArgs.counterMax
                NameValueArgs.dsAfter
            end
            if isempty(NameValueArgs.counterMax)
                NameValueArgs = rmfield(NameValueArgs,'counterMax');
            end
            if isempty(NameValueArgs.dsAfter)
                NameValueArgs = rmfield(NameValueArgs,'dsAfter');
            end
            if isempty(obj.allEvents)
                obj.allEvents = {stepSize.StepSizeSingleEvent(eventCondition, neededParams, dsMin, dsMax, NameValueArgs)};
            else
                len = numel(obj.allEvents);
                tmp = cell(len+1,1);
                for k = 1:len
                    tmp{k} = obj.allEvents{k};
                end
                tmp{end} = stepSize.StepSizeSingleEvent(eventCondition, neededParams, dsMin, dsMax, NameValueArgs);
                obj.allEvents = tmp;
            end
        end

        function output = getActiveEvents(obj)
            output = obj.activeEvents;
        end
        
        function output = getChanged(obj)
            output = obj.changed;
        end

        function output = getNeededVariables(obj)
            if isempty(obj.allParams)
                list = {};
                for k = 1:numel(obj.allEvents)
                    params = obj.allEvents{k}.get('neededParameters');
                    list = [list, params];
                end
                obj.allParams = unique(list);
            end
            output = obj.allParams;
        end

        function obj = setInitial(obj,Initial)
            obj.Initial = Initial;
        end

        function obj = setCurrent(obj,dsMax,dsMin)
            obj.Current.dsMax = dsMax;
            obj.Current.dsMin = dsMin;
        end
    end
    methods (Access = private)
        function newStepsizes = setStepsize(obj,isActive,dsCurrent)
            dominatingEvent = NaN;
            ds = inf;
            newStepsizes.dsMax = [];
            newStepsizes.dsMin = [];
            obj.changed = false;
            % iterate through events
            for k = 1:numel(obj.allEvents)
                % if any is active and its stepsize is lower than the
                % previous saved value and it wasn't active in the latest step
                % save its value
                if isActive(k) && ~obj.allEvents{k}.get('lastActive') && obj.allEvents{k}.get('dsMax') < obj.Current.dsMax
                    obj.changed = true;
                    ds = obj.allEvents{k}.get('dsMax');
                    obj.allEvents{k}.set('lastActive',true);

                    % set min and max stepsize
                    newStepsizes.dsMax = obj.allEvents{k}.get('dsMax');
                    newStepsizes.dsMin = obj.allEvents{k}.get('dsMin');

                    % save event number
                    dominatingEvent = k;
                end
            end
            %
            % check whether found stepsize is bigger than currently used one
            %
            if ds >= dsCurrent
                newStepsizes.ds = dsCurrent;
            else
                newStepsizes.ds = ds;
            end
            %
            % increase event counter
            %
            if ~isnan(dominatingEvent)
                obj.allEvents{dominatingEvent}.increaseCounter();
            end
        end
    end
end