%% path continuation - step_size.event
%
%   Institute of Dynamics and Vibration Research
%   Leibniz University Hannover
%   21.01.2022 - Tido Kubatschek
%
%
function [ds,Path] = event(Path,Opt)
    % Benoetigt:
    % Opt-Eintraege, ob es a
    if eval(Opt.stepsize_event) && Path.event_flag
        ds = Opt.event_stepsize;
        Path.event_flag = 0;
    end
end