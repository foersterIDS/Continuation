%% path continuation - closed_curve
%
%   Institute of Dynamics and Vibration Research
%   Leibniz University Hannover
%   03.11.2020 - Tido Kubatschek
%   09.11.2020 - Alwin Förster
%
function [is_closed] = closed_curve(Opt, var_all,l_all, s_all, ds)
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

        flag = 0;
        if numel(l_all)>n
            dist_x = sqrt(sum((x_all(:,nu:numel(l_all) - 2)-x_all(:,end)).^2,1));
            k_flags = find(dist_x <= eps_dist);
            if numel(k_flags)>0
                flag = 1;
                k_flags = k_flags + nu - 1;
            end
        end

        %% check wether the found point is crossed under the same angle
        if flag == 1
            for i=1:numel(k_flags)
                %% found point (f) and the one before (fb) and after (fa)
                vec_f = x_all(:,k_flags(i)) - x_all(:,k_flags(i)-1);
                s_f = s_all(k_flags(i));
                s_fb = s_all(k_flags(i) - 1);
                s_fa = s_all(k_flags(i) + 1);
                point_f = x_all(:,k_flags(i));
                point_fb = x_all(:,k_flags(i)-1);
                point_fa = x_all(:,k_flags(i)+1);

                %% determine parabola for those three points in regards to s
                p_f = polyfitn([s_fb, s_f, s_fa], [point_fb, point_f, point_fa],2);
                %% current point (c) and one before (cb) and two steps before (cbb)
                vec_c = x_all(:,end) - x_all(:,end - 1);
                s_c = s_all(end );
                s_cb = s_all(end - 1);
                s_cbb = s_all(end - 2);
                point_c = x_all(:,end);
                point_cb = x_all(:,end - 1);
                point_cbb = x_all(:,end - 2);

                %% parabola for those three points in regards to s
                p_c = polyfitn([s_cbb, s_cb, s_c], [point_cbb, point_cb, point_c],2);

                % shifting parabolas to origin and therefore just need
                % to compare first parameters:

                if abs(p_c(1)) >= abs(p_f(1)) %always devide by bigger one
                    comp = (p_c(1)-p_f(1)) / p_c(1);
                else
                    comp = (p_c(1)-p_f(1)) / p_f(1);
                end

                if abs(abs(comp) - 1) <= 1e-2
                    Opt.closed_counter = Opt.closed_counter - 1;
                end

%                     %% angle
%                     angle = vector_angle(vec_f, vec_c);
%                     %% direction
%                     r_f = vec_f / norm(vec_f);
%                     r_c = vec_c / norm(vec_c);
%                     dir = norm(r_f - r_c);
%                     %% check
%                     if angle <= eps_ang && dir <= eps_dir && comp > 0.9 % check angle and direction
%                         closed_counter = closed_counter + 1;
%                     end
            end
        end

        if Opt.closed_counter == 0
            is_closed = 1;
        end
    end
end

