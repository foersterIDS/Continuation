%% path continuation - step_size.event_adjustment
%
%  See the <a href="matlab:open('..\doc\html\continuation.html')">documentation</a> or see <a href="matlab:doc('step_size.control')">other stepsize adaption methods</a>.
%  See <a href="matlab:doc('StepSizeControlEvent')">StepSizeControlEvent class</a> for information on how to implement events
%
%   Institute of Dynamics and Vibration Research
%   Leibniz University Hannover
%   08.11.2022 - Tido Kubatschek
%
function [ds_new,event_obj,changed,Opt] = event_adjustment(ds,event_obj,Opt,Initial,variables)
    %
    % check if event_object has already been defined
    if isempty(event_obj)
        %
        % create event_obj
        %
        event_obj = step_size.StepSizeEventManager;
        %
        % get information by user
        %
        len = length(Opt.event_user_input);
        %
        % check len
        %
        if len == 0
            error('step_size_event has been activated, but there are no entries in user input!');
        end
        %
        for k = 1:len
            %
            % get inputs
            %
            condition = Opt.event_user_input{k}.condition;
            needed_parameters = Opt.event_user_input{k}.needed_parameters;
            ds_min = Opt.event_user_input{k}.ds_min;
            ds_max = Opt.event_user_input{k}.ds_max;

            if isfield(Opt.event_user_input{k},'counter_max')
                var_input = {'counter',Opt.event_user_input{k}.counter_max};
            else
                var_input = [];
            end
            %
            % create event
            %
            event_obj.add_event(condition,needed_parameters,ds_min,ds_max,var_input);
            %
        end
        %
        % set initial ds_min and ds_max
        %
        event_obj.set_initial(Initial);
        %
    end
    %
    % get list of needed variables
    needed_parameters = event_obj.getNeededVariables;
    %
    % set current ds_max and ds_min
    %
    event_obj.set_current(Opt.ds_max,Opt.ds_min); % ist hier Opt richtig?
    %
    % fill value struct
    num_of_vars = length(needed_parameters);
    for k = 1:num_of_vars
        eval(['value_struct.',needed_parameters{k},'=variables.',needed_parameters{k},';']);
    end
    %
    current_stepsize = ds;
    new_stepsizes = event_obj.check_events_and_adapt_stepsize(current_stepsize,value_struct);
    changed = event_obj.getChanged;
    %
    % adapt ds_max, ds_min and ds
    %
    if ~isempty(new_stepsizes.ds_max)
        event_obj.set_current(new_stepsizes.ds_max,new_stepsizes.ds_min);
        Opt.ds_max = new_stepsizes.ds_max;
        Opt.ds_min = new_stepsizes.ds_min;
    end
    if changed
        ds_new = new_stepsizes.ds;
    else
        ds_new = ds;
    end
    %
end