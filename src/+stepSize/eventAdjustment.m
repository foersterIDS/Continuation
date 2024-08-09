%% path continuation - stepSize.eventAdjustment
%
%  See the <a href="matlab:open('..\doc\html\continuation.html')">documentation</a> or see <a href="matlab:doc('stepSize.control')">other stepsize adaption methods</a>.
%  See <a href="matlab:doc('StepSizeControlEvent')">StepSizeControlEvent class</a> for information on how to implement events
%
%   Institute of Dynamics and Vibration Research
%   Leibniz University Hannover
%   08.11.2022 - Tido Kubatschek
%
function [dsNew,eventObj,changed] = eventAdjustment(ds,eventObj,oih,variables)
    %
    % check if eventObject has already been defined
    if isempty(eventObj)
        %
        % create eventObj
        %
        eventObj = stepSize.StepSizeEventManager;
        %
        % get information by user
        %
        len = length(oih.opt.eventUserInput);
        %
        % check len
        %
        if len == 0
            error('stepSizeEvent has been activated, but there are no entries in user input!');
        end
        %
        for k = 1:len
            %
            % get inputs
            %
            condition = oih.opt.eventUserInput{k}.condition;
            neededParameters = oih.opt.eventUserInput{k}.neededParameters;
            dsMin = oih.opt.eventUserInput{k}.dsMin;
            dsMax = oih.opt.eventUserInput{k}.dsMax;
            
            if isfield(oih.opt.eventUserInput{k},'counterMax')
                counterMax = oih.opt.eventUserInput{k}.counterMax;
            else
                counterMax = [];
            end

            if isfield(oih.opt.eventUserInput{k},'dsAfter')
                dsAfter = oih.opt.eventUserInput{k}.dsAfter;
            else
                dsAfter = [];
            end
            %
            % create event
            %
            eventObj.addEvent(condition,neededParameters,dsMin,dsMax,"counterMax",counterMax,"dsAfter",dsAfter);
            %
        end
        %
        % set initial dsMin and dsMax
        %
        eventObj.setInitial(oih.initial);
        %
    end
    %
    % get list of needed variables
    neededParameters = eventObj.getNeededVariables;
    %
    % set current dsMax and dsMin
    %
    eventObj.setCurrent(oih.opt.dsMax,oih.opt.dsMin); % ist hier Opt richtig?
    %
    % fill value struct
    numOfVars = length(neededParameters);
    for k = 1:numOfVars
        eval(['valueStruct.',neededParameters{k},'=variables.',neededParameters{k},';']);
    end
    %
    currentStepsize = ds;
    newStepsizes = eventObj.checkEventsAndAdaptStepsize(currentStepsize,valueStruct);
    changed = eventObj.getChanged;
    %
    % adapt dsMax, dsMin and ds
    %
    if ~isempty(newStepsizes.dsMax)
        eventObj.setCurrent(newStepsizes.dsMax,newStepsizes.dsMin);
        oih.opt.dsMax = newStepsizes.dsMax;
        oih.opt.dsMin = newStepsizes.dsMin;
    end
    if changed
        dsNew = newStepsizes.ds;
    else
        dsNew = ds;
    end
    %
end