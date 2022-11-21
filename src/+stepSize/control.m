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
function [dsn,Counter,event,Opt] = control(ds,Counter,Solver,Do,Plus,Path,Jacobian,Opt,Info,event,Initial)
    if ~Do.stepback
        if ~Do.deflate
            if Counter.error == 0
                if Opt.stepSizeEvent                 
                    % fill variables struct
                    % in the first run, fill with all possible variables
                    if isempty(event.eventObj)
                        event.variables.dsCurrent = ds;
                        event.variables.dsMinCurrent = Opt.dsMin;
                        event.variables.dsMaxCurrent = Opt.dsMax;
                        event.variables.Counter = Counter;
                        event.variables.Solver = Solver;
                        event.variables.Path = Path;
                        event.variables.Jacobian = Jacobian;
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
                                    event.variables.dsMinCurrent = Opt.dsMin;
                                case 'dsMaxCurrent'
                                    event.variables.dsMaxCurrent = Opt.dsMax;
                                otherwise
                            event.variables.(event.namesOfNeededVariables{k}) = eval(event.namesOfNeededVariables{k});
                            end
                        end
                    end

                    % adjust stepsize
                    [dsn,event.eventObj,changed,Opt] = stepSize.eventAdjustment(ds,event.eventObj,Opt,Initial,event.variables);
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
                    if Opt.stepSizeControl.angleChange
                        %
                        % Check if there are enough solution points to use
                        % method. Otherwise use iterations.
                        %
                        if length(Path.lAll) > 3
                            xi = stepSize.angleChange(Solver,Path,Opt);
                        else
                            xi = stepSize.iterationsPolynomial(Solver,Opt);
                        end
                    %
                    % angle custom
                    %
                    elseif Opt.stepSizeControl.angleCustom
                        %
                        % Check if there are enough solution points to use
                        % method. Otherwise use iterations.
                        %
                        if length(Path.lAll) > 2
                            xi = stepSize.angleCustom(Solver,Path,Opt);
                        else
                            xi = stepSize.iterationsPolynomial(Solver,Opt);
                        end
                    %
                    % contraction
                    %
                    elseif Opt.stepSizeControl.contraction
                        xi = stepSize.contraction(Solver,Opt);
                    %
                    % error
                    %
                    elseif Opt.stepSizeControl.error
                        xi = stepSize.error(Solver,Path,Opt);
                    %
                    % errorAlt
                    %
                    elseif Opt.stepSizeControl.errorAlt
                        xi = stepSize.errorAlt(Solver,Path,Opt);
                    %
                    % fayezioghani
                    %
                    elseif Opt.stepSizeControl.fayezioghani
                        %
                        % Check if there are enough solution points to use
                        % method. Otherwise use iterations.
                        %
                        if length(Path.lAll) > 2
                            xi = stepSize.fayezioghani(ds,Solver,Path,Jacobian.last,Opt);
                        else
                            xi = stepSize.iterationsPolynomial(Solver,Opt);
                        end
                    %
                    % fixed step size
                    %
                    elseif Opt.stepSizeControl.fix
                        xi = 1;
                    % 
                    % iterations based - exponential
                    %
                    elseif Opt.stepSizeControl.iterationsExponential
                        xi = stepSize.iterationsExponential(Solver,Opt);
                    % 
                    % iterations based - polynomial
                    %
                    elseif Opt.stepSizeControl.iterationsPolynomial
                        xi = stepSize.iterationsPolynomial(Solver,Opt);
                    % 
                    % multiplicative method
                    %
                    elseif Opt.stepSizeControl.multiplicative
                        xi = stepSize.multiplicative(Solver,Path,Opt);
                    % 
                    % multiplicativeAlt method
                    %
                    elseif Opt.stepSizeControl.multiplicativeAlt
                        xi = stepSize.multiplicativeAlt(Solver,Path,Opt);
                    %
                    % pid control - custom
                    %
                    elseif Opt.stepSizeControl.pidCustom
                        %
                        % Check if there are enough solution points to use
                        % method. Otherwise use iterations.
                        %
                        if length(Path.lAll) > 4
                            xi = stepSize.pidCustom(Path,Opt);
                        else
                            xi = stepSize.iterationsPolynomial(Solver,Opt);
                        end
                    %
                    % pid control - valli
                    %
                    elseif Opt.stepSizeControl.pidValli
                        %
                        % Check if there are enough solution points to use
                        % method. Otherwise use iterations.
                        %
                        if length(Path.lAll) > 4
                            xi = stepSize.pidValli(Path,Opt);
                        else
                            xi = stepSize.iterationsPolynomial(Solver,Opt);
                        end
                    %
                    % szyszkowski
                    %
                    elseif Opt.stepSizeControl.szyszkowski
                        if length(Path.lAll) > 3
                            xi = stepSize.szyszkowski(Solver,Path,Opt);
                        else
                            xi = stepSize.iterationsPolynomial(Solver,Opt);
                        end
                    %
                    % yoon
                    %
                    elseif Opt.stepSizeControl.yoon
                        if length(Path.lAll) > 3
                            xi = stepSize.yoon(Solver,Path,Opt);
                        else
                            xi = stepSize.iterationsPolynomial(Solver,Opt);
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
                    %% Limit to max./min. step size, also limit to Opt.maxStepSizeChange*ds / ds/Opt.maxStepSizeChange*ds:
                    dsn = min([norm(Opt.dsMax),Opt.maxStepSizeChange*ds,dsn]);
                    dsn = max([Opt.dsMin,ds/Opt.maxStepSizeChange,dsn]);
                end
            else
                if Opt.stepSizeControl.fix
                    dsn = Info.ds0;
                else
                    dsn = ds/2;
                    %% Limit to max./min. step size, also limit to Opt.maxStepSizeChange*ds / ds/Opt.maxStepSizeChange*ds:
                    dsn = min([norm(Opt.dsMax),Opt.maxStepSizeChange*ds,dsn]);
                    dsn = max([Opt.dsMin,ds/Opt.maxStepSizeChange,dsn]);
                end
            end
        else
            dsn = ds;
        end
    else
        if Opt.stepSizeControl.fix
            dsn = Info.ds0;
        else
            xe = [Path.varAll(:,end);Path.lAll(end)];
            dsn = norm(Plus.x-xe)/2;
        end
    end
end