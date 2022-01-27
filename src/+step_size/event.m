%% path continuation - step_size.event
%
%   Institute of Dynamics and Vibration Research
%   Leibniz University Hannover
%   21.01.2022 - Tido Kubatschek
%
%   Inputs:
%       ds          -- latest stepsize
%       Path        -- path information
%       Counter     -- contains event_counter
%       Opt         -- contains user inputs and preferences such as
%                      ds_event
%       event_out   -- true --> event condition was met in the last step 
%                               and stepsize has been adapted. The step
%                               size can only be adjusted again if the
%                               condition is not met at least once.                                                  
%                      false --> event condition was not met in the last 
%                                step, so it can be adapted again.
%                        
%   Outputs:
%       ds          -- adapted stepsize / latest stepsize
%       Counter     -- contains event_counter
%       event_out   -- ...
%       changed     -- tells if the stepsize has been changed
%
function [ds,Counter,event_out,changed] = event(ds,Path,Counter,Opt,event_out)
    %
    % To-Do: identify relevant variables for event checking
    %        --> Path, ds, Jacobian?, fun?
    %
    %
    % Check if stepsize is still allowed to be changed
    %
    changed = false;
    if Opt.event_condition(ds,Path) && Counter.event < Opt.event_counter
        %
        % check if condition was met in the last step
        %
        if ~event_out
            %
            % stepsize is allowed to be changed now
            ds = Opt.ds_event;
            Counter.event = Counter.event + 1;
            changed = true;
            %
        end
            %
            % stepsize isnt allowed to be changed in the next step (since
            % event occured in the last / this step)
            %
            event_out = true;
    else
        event_out = false;
    end
end