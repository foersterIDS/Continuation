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
        flag = 0;
        if numel(l_all)>n
            dist_x = sqrt(sum((x_all(:,nu:numel(l_all)-2)-x_all(:,end)).^2,1));
            k_flags = find(dist_x<=eps_dist);
            if numel(k_flags)>0
                flag = 1;
                k_flags = k_flags+nu-1;
            end
        end
        
        %% check wether the found point is crossed under the same angle
        if flag == 1
            for i=1:numel(k_flags)
                %% found point
                vec_f = x_all(:,k_flags(i)) - x_all(:,k_flags(i)-1);
                %% current point
                vec_c = x_all(:,end) - x_all(:,end-1);
                %% angle
                angle = vector_angle(vec_f, vec_c);
                %% normalized scalar prod
                normprod = dot(vec_f, vec_f) / (norm(vec_f) * norm(vec_c));
                %% direction
                r_f = vec_f / norm(vec_f);
                r_c = vec_c / norm(vec_c);
                dir = norm(r_f - r_c);
                %% ceck
                if angle <= eps_ang % check angle
                    is_closed = 1;
                    break;
%                 elseif abs(1-normprod) <= 1e-6 % check normalized scalar product
%                     is_closed = 1;
%                     break;
                elseif dir <= eps_dir % check direction
                    is_closed = 1;
                    break;
                end
            end
        end
    end
end

