%% path continuation - step_size_control_angle
%
%   Institute of Dynamics and Vibration Research
%   Leibniz University Hannover
%   07.10.2020 - Tido Kubatschek
%   08.10.2020 - Alwin Förster
%
function [dsn] = step_size_control_angle(ds,ds0,error_counter,solver_output,do_deflate,vars,ls,Opt)
    %% calc angular change
    if length(ls) > 2
%         vh = length(vars(:,1));
%         theta = 0;
%         
%         % calculate greatest angle
%         for k = 1:vh
%             vr = vars(k,end-2:end);
%             lr = ls(end-2:end);
%             vec = [lr;vr];
%             v1 = vec(:,end) - vec(:,end-1);
%             v2 = vec(:,end-1) - vec(:,end-2);
%             
%             thetn = real(acos(max(min(dot(v1,v2)/(norm(v1)*norm(v2)),1),-1)));
% 
%             if thetn > theta
%                 theta = thetn;
%             end
%         end
        % Vorschlag für Winkelberechnung:
        vec = [vars;ls];
        v1 = vec(:,end) - vec(:,end-1);
        v2 = vec(:,end-1) - vec(:,end-2);
        theta = acos(dot(v1,v2)/(norm(v1)*norm(v2)));
        
        % calculate ratio
        ch = theta / Opt.step_size_angle;
        
        if ch < 1 % if ratio < 1, set ch to 1
            ch = 1;
        elseif ch > 1.5 % high ratios are punished stronger (but ch_max = 5)
            ch = 5;
        end
    else
        ch = 1;
    end
    % calculate step size
    dsn = ds*Opt.n_iter_opt/(solver_output.iterations*ch);
    dsn = max(ds/2,dsn);
    dsn = min(ds*2,dsn);
end