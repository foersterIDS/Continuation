%% path continuation - step_size.angle
%
%   Institute of Dynamics and Vibration Research
%   Leibniz University Hannover
%   07.10.2020 - Tido Kubatschek
%   08.10.2020 - Alwin Förster
%
function [dsn] = angle(ds,ds0,Counter,solver_output,Do,Path,Opt)
    %% calc angular change
    if length(Path.l_all) > 2
        vh = length(Path.var_all(:,1));
        angle_max = 0;
        
        % calculate greatest angle
        for k = 1:vh
            vr = Path.var_all(k,end-2:end);
            lr = Path.l_all(end-2:end);
            vec = [lr;vr];
            v1 = vec(:,end) - vec(:,end-1);
            v2 = vec(:,end-1) - vec(:,end-2);
            angle_n = aux.vector_angle(v1,v2);

            if angle_n > angle_max
                angle_max = angle_n;
            end
        end
        
        angle = angle_max;
        
        % calculate ratio
        ch = angle / Opt.step_size_angle;
        
        if ch < 1 % if ratio < 1, set ch to 1
            ch = 1;
        elseif ch > 1.5 && ch < 5 % high ratios are punished stronger (but ch_max = 5)
            ch = 5;
        elseif ch > 5
            fprintf('--> Point of high curvature!\n');
            ch = 10;
        end
    else
        ch = 1;
    end
    
    % calculate step size
    dsn = ds*sqrt(Opt.n_iter_opt/(solver_output.iterations*ch));
end