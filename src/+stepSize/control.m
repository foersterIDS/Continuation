%% path continuation - stepSize.control
%  Adjusts stepsize by choosing an adaption method. 
%
%
%
%  See the <a href="matlab:open('..\doc\html\continuation.html')">documentation</a> or
%  see the control methods:
%  -- <a href="matlab:doc('stepSize.angleChange')">stepSize.angleChange</a>
%  -- <a href="matlab:doc('stepSize.angleCustom')">stepSize.angleCustom</a>
%  -- <a href="matlab:doc('stepSize.contraction')">stepSize.contraction</a>
%  -- <a href="matlab:doc('stepSize.error')">stepSize.error</a>
%  -- <a href="matlab:doc('stepSize.eventAdjustment')">stepSize.eventAdjustment</a>
%  -- <a href="matlab:doc('stepSize.fayezioghani')">stepSize.fayezioghani</a>
%  -- <a href="matlab:doc('stepSize.iterationsExponential')">stepSize.iterationsExponential</a>
%  -- <a href="matlab:doc('stepSize.iterationsPolynomial')">stepSize.iterationsPolynomial</a>
%  -- <a href="matlab:doc('stepSize.multiplicative')">stepSize.multiplicative</a>
%  -- <a href="matlab:doc('stepSize.pidCustom')">stepSize.pidCustom</a>
%  -- <a href="matlab:doc('stepSize.pidValli')">stepSize.pidValli</a>
%  -- <a href="matlab:doc('stepSize.szyszkowski')">stepSize.szyszkowski</a>
%  -- <a href="matlab:doc('stepSize.yoon')">stepSize.yoon</a>
%
%
%   Institute of Dynamics and Vibration Research
%   Leibniz University Hannover
%   08.05.2020 - Alwin Förster
%   03.06.2020 - Niklas Marhenke
%   21.10.2020 - Tido Kubatschek
%
function [dsn,event] = control(ds,oih,event)
    if ~oih.do.stepback
        if ~oih.do.deflate
            if oih.counter.error == 0
                if oih.opt.stepSizeEvent                 
                    % fill variables struct
                    % in the first run, fill with all possible variables
                    if isempty(event.eventObj)
                        event.variables.dsCurrent = ds;
                        event.variables.dsMinCurrent = oih.opt.dsMin;
                        event.variables.dsMaxCurrent = oih.opt.dsMax;
                        event.variables.Counter = oih.counter;
                        event.variables.Solver = oih.solver;
                        event.variables.Path = oih.path;
                    else
                        % after creating eventObj, read out names of
                        % needed variables with method and save them in
                        % field of event
                        if ~isfield(event,'namesOfNeededVariables')
                            % get names of needed variables
                            event.namesOfNeededVariables = event.eventObj.getNeededVariables;
                            % empty event.variables
                            event.variables = [];
                        end
                        
                        % loop through names of needed variables and save
                        % them on event.variables
                        for k = 1:length(event.namesOfNeededVariables)
                            switch event.namesOfNeededVariables{k}
                                case 'dsCurrent'
                                    event.variables.dsCurrent = ds;
                                case 'dsMinCurrent'
                                    event.variables.dsMinCurrent = oih.opt.dsMin;
                                case 'dsMaxCurrent'
                                    event.variables.dsMaxCurrent = oih.opt.dsMax;
                                otherwise
                            event.variables.(event.namesOfNeededVariables{k}) = eval(event.namesOfNeededVariables{k});
                            end
                        end
                    end

                    % adjust stepsize
                    [dsn,event.eventObj,changed,Opt] = stepSize.eventAdjustment(ds,event.eventObj,oih,event.variables);
                else
                    changed = false;
                end
                if ~changed
                    %
                    %% Determine step size control method
                    %
                    %
                    % angle change based
                    %
                    if oih.opt.stepSizeControl.angleChange
                        %
                        % Check if there are enough solution points to use
                        % method. Otherwise use iterations.
                        %
                        if oih.path.nAll > 3
                            xi = stepSize.angleChange(oih);
                        else
                            xi = stepSize.iterationsPolynomial(oih);
                        end
                    %
                    % angle custom
                    %
                    elseif oih.opt.stepSizeControl.angleCustom
                        %
                        % Check if there are enough solution points to use
                        % method. Otherwise use iterations.
                        %
                        if oih.path.nAll > 2
                            xi = stepSize.angleCustom(oih);
                        else
                            xi = stepSize.iterationsPolynomial(oih);
                        end
                    %
                    % contraction
                    %
                    elseif oih.opt.stepSizeControl.contraction
                        xi = stepSize.contraction(oih);
                    %
                    % error
                    %
                    elseif oih.opt.stepSizeControl.error
                        xi = stepSize.error(oih);
                    %
                    % errorAlt
                    %
                    elseif oih.opt.stepSizeControl.errorAlt
                        xi = stepSize.errorAlt(oih);
                    %
                    % fayezioghani
                    %
                    elseif oih.opt.stepSizeControl.fayezioghani
                        %
                        % Check if there are enough solution points to use
                        % method. Otherwise use iterations.
                        %
                        if oih.path.nAll > 2
                            xi = stepSize.fayezioghani(ds,oih);
                        else
                            xi = stepSize.iterationsPolynomial(oih);
                        end
                    %
                    % fixed step size
                    %
                    elseif oih.opt.stepSizeControl.fix
                        xi = 1;
                    % 
                    % iterations based - exponential
                    %
                    elseif oih.opt.stepSizeControl.iterationsExponential
                        xi = stepSize.iterationsExponential(oih);
                    % 
                    % iterations based - polynomial
                    %
                    elseif oih.opt.stepSizeControl.iterationsPolynomial
                        xi = stepSize.iterationsPolynomial(oih);
                    % 
                    % multiplicative method
                    %
                    elseif oih.opt.stepSizeControl.multiplicative
                        xi = stepSize.multiplicative(oih);
                    % 
                    % multiplicativeAlt method
                    %
                    elseif oih.opt.stepSizeControl.multiplicativeAlt
                        xi = stepSize.multiplicativeAlt(oih);
                    %
                    % pid control - custom
                    %
                    elseif oih.opt.stepSizeControl.pidCustom
                        %
                        % Check if there are enough solution points to use
                        % method. Otherwise use iterations.
                        %
                        if oih.path.nAll > 4
                            xi = stepSize.pidCustom(oih);
                        else
                            xi = stepSize.iterationsPolynomial(oih);
                        end
                    %
                    % pid control - valli
                    %
                    elseif oih.opt.stepSizeControl.pidValli
                        %
                        % Check if there are enough solution points to use
                        % method. Otherwise use iterations.
                        %
                        if oih.path.nAll > 4
                            xi = stepSize.pidValli(oih);
                        else
                            xi = stepSize.iterationsPolynomial(oih);
                        end
                    %
                    % szyszkowski
                    %
                    elseif oih.opt.stepSizeControl.szyszkowski
                        if oih.path.nAll > 3
                            xi = stepSize.szyszkowski(oih);
                        else
                            xi = stepSize.iterationsPolynomial(oih);
                        end
                    %
                    % yoon
                    %
                    elseif oih.opt.stepSizeControl.yoon
                        if oih.path.nAll > 3
                            xi = stepSize.yoon(oih);
                        else
                            xi = stepSize.iterationsPolynomial(oih);
                        end
                    %
                    % wrong method
                    %
                    else
                        error('Invalid settings for step size control!');
                    end
                    %
                    %% adapt step size
                    dsn = xi * ds;
                    %% Limit to max./min. step size, also limit to oih.opt.maxStepSizeChange*ds / ds/oih.opt.maxStepSizeChange*ds:
                    dsn = min([norm(oih.opt.dsMax),oih.opt.maxStepSizeChange*ds,dsn]);
                    dsn = max([oih.opt.dsMin,ds/oih.opt.maxStepSizeChange,dsn]);
                end
            else
                if oih.opt.stepSizeControl.fix
                    dsn = oih.info.ds0;
                else
                    dsn = ds/2;
                    %% Limit to max./min. step size, also limit to oih.opt.maxStepSizeChange*ds / ds/oih.opt.maxStepSizeChange*ds:
                    dsn = min([norm(oih.opt.dsMax),oih.opt.maxStepSizeChange*ds,dsn]);
                    dsn = max([oih.opt.dsMin,ds/oih.opt.maxStepSizeChange,dsn]);
                end
            end
        else
            dsn = ds;
        end
    else
        if oih.opt.stepSizeControl.fix
            dsn = oih.info.ds0;
        else
            xe = [oih.path.varAll(:,end);oih.path.lAll(end)];
            dsn = norm(oih.path.xPlus-xe)/2;
        end
    end
end