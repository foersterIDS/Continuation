%% path continuation - closed_curve
%
%   Institute of Dynamics and Vibration Research
%   Leibniz University Hannover
%   03.11.2020 - Tido Kubatschek
%
function [is_closed] = closed_curve(var_all,l_all, ds)
    is_closed = 0;
    %letze Punkt nicht betrachten siehe 2 fache Schrittweite
    % ggf fragen, ob trotzdem weitergemacht werden soll
    n = 5;
    if length(l_all) > n+2
        eps_dist = ds; %Should be prop to ds to avoid errors
        eps_ang = 10*(2*pi / 360);

        x_all = [var_all; l_all];

        %% search in x_all if last solution has already been found
        flag = 0;
        for k = n:length(l_all)-1-n
            dist_x = norm(x_all(:,k) - x_all(:,end));
            if dist_x <= eps_dist
                flag = 1;
                k_flag = k;
            end
        end

        %% check wether the found point is crossed under the same angle
        if flag == 1
            %% found point
            vec_f = x_all(:,k_flag) - x_all(:,k_flag-1);
            %% current point
            vec_c = x_all(:,end) - x_all(:,end-1);

            %% angle
            angle = vector_angle(vec_f, vec_c);
            
            %% normalized scalar prod
            normprod = dot(vec_f, vec_f) / (norm(vec_f) * norm(vec_c));

            %% ceck angle
            if angle <= eps_ang
                is_closed = 1;
            end
            
            %% check normalized scalar product
            if abs(1-normprod) <= 1e-6
                is_closed = 1;
            end
        end
    end
end

