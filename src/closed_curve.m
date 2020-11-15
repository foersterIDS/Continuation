%% path continuation - closed_curve
%
%   Institute of Dynamics and Vibration Research
%   Leibniz University Hannover
%   03.11.2020 - Tido Kubatschek
%   09.11.2020 - Alwin Förster
%
function [is_closed] = closed_curve(var_all,l_all, ds)
    is_closed = 0;
    %letze Punkt nicht betrachten siehe 2 fache Schrittweite
    % ggf fragen, ob trotzdem weitergemacht werden soll
    n = 5;
    nu = 2;
    if length(l_all) > n+nu
        x_all = [var_all; l_all];
        eps_dist = norm(x_all(:,end)-x_all(:,end-1));
        eps_ang = 1*(2*pi / 360);
        eps_dir = 1e-1;
        %% search in x_all if last solution has already been found
        closed_counter = 0;
        
        % how many of the latest points must be found a similar one for
        closed_counter_max = 3;
        
        for kk=0:closed_counter_max-1
            flag = 0;
            if numel(l_all)>n
                dist_x = sqrt(sum((x_all(:,nu:numel(l_all)- kk - 2)-x_all(:,end - kk)).^2,1));
                k_flags = find(dist_x <= eps_dist);
                if numel(k_flags)>0
                    flag = 1;
                    k_flags = k_flags + nu - 1;
                end
            end

            %% check wether the found point is crossed under the same angle
            if flag == 1
                for i=1:numel(k_flags)
                    %% found point
                    vec_f = x_all(:,k_flags(i)) - x_all(:,k_flags(i)-1);
                    %% current point
                    vec_c = x_all(:,end - kk) - x_all(:,end - 1 -kk);
                    %% angle
                    angle = vector_angle(vec_f, vec_c);
                    %% direction
                    r_f = vec_f / norm(vec_f);
                    r_c = vec_c / norm(vec_c);
                    dir = norm(r_f - r_c);
                    %% check
                    if angle <= eps_ang && dir <= eps_dir % check angle and direction
                        closed_counter = closed_counter + 1;
                        break;
                    end
                end
            end
        end
        if closed_counter == closed_counter_max
            is_closed = 1;
        end
    end
end

