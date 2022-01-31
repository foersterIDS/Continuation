%% path continuation - aux.closed_curve
%
%   Institute of Dynamics and Vibration Research
%   Leibniz University Hannover
%   03.11.2020 - Tido Kubatschek
%   09.11.2020 - Alwin Förster
%
function [is_closed, Opt, Counter] = closed_curve(Opt,Path,ds,Counter)
    is_closed = 0;
    n = 5;
    np = 4; % polynomial order (np = 2*k with k in N)
    eps_pol = 2.5e-2;
    nu = 1+np/2;
    ns = numel(Path.l_all);
    if ns > n+nu
        x_all = [Path.var_all; Path.l_all];
        eps_dist = min([norm(x_all(:,end)-x_all(:,end-1)),2*ds,Opt.ds_max]);
        eps_ang = 0.5*(2*pi / 360);
        eps_dir = 1e-3;
        
        %% exclude points before current points which are too close
        %
        ignored_dist = 2*norm(x_all(:,end)-x_all(:,end-1));
        distance = sqrt(sum((x_all(:,1:(ns - 1)) - x_all(:,end)).^2,1));
        last_ind = find(distance > ignored_dist);
        if isempty(last_ind)
            flag = 0;
        else
            ignored = ns - last_ind(end);
            %
            flag = 0;
            if ns>n
                dist_x = distance(nu:(ns - ignored));
                k_flags = find(dist_x <= eps_dist);
                if numel(k_flags)>0
                    flag = 1;
                    k_flags = k_flags + nu - 1;
                end
            end
        end       
        %% check wether the found point is crossed under the same angle
        count_flag = 0;
        if flag == 1
            for i=1:numel(k_flags)
                %% polynomials
                % centered polynomial to matching point
                s0_f = Path.s_all(k_flags(i));
                p_f = poly.fitn(Path.s_all(k_flags(i)+(-np/2:np/2))-s0_f,x_all(:,k_flags(i)+(-np/2:np/2)),np);
                % arclength correction
                dsc = dist_x(k_flags(i)-nu+1);
                % centered polynomials to last point
                s0_c = Path.s_all(end);
                p_c_p = poly.fitn(Path.s_all(end+(-np:0))-s0_c+dsc,x_all(:,end+(-np:0)),np);
                p_c_m = poly.fitn(Path.s_all(end+(-np:0))-s0_c-dsc,x_all(:,end+(-np:0)),np);
                % check polyinomials
                if min([norm(p_f-p_c_p),norm(p_f-p_c_m)]/norm(p_f)) <= eps_pol
                    count_flag = 1;
                    break;
                end
                %% angle
                vec_f = x_all(:,k_flags(i)) - x_all(:,k_flags(i)-1);
                vec_c = x_all(:,end) - x_all(:,end - 1);
                angle = aux.vector_angle(vec_f, vec_c);
                %% direction
                r_f = vec_f / norm(vec_f);
                r_c = vec_c / norm(vec_c);
                dir = norm(r_f - r_c);
                %% check
                if angle <= eps_ang && dir <= eps_dir % check angle and direction
                    count_flag = 1;
                end
            end
        end
        if count_flag
            Counter.closed_curve = Counter.closed_curve + 1;
        elseif Counter.closed_curve > 0
            Counter.closed_curve = Counter.closed_curve - 1;
        end
        if Counter.closed_curve == Opt.max_closed_counter
            is_closed = 1;
        end
    end
end

