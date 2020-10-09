%% path continuation - live_plot
%
%   Institute of Dynamics and Vibration Research
%   Leibniz University Hannover
%   30.09.2020 - Tido Kubatschek
%
function [pl_info] = live_plot(Opt, nv, l_start, l_end, l_all, var_all, pl_info, bif_flag, bif)
    l_lu = [min([l_start,l_end]),max([l_start,l_end])];
    l_max = [min(l_all),max(l_all)];
    dl0 = abs(l_end-l_start)*0.05;
    if length(l_all) == 1
        fig = figure('units', 'normalized', 'position', [0.2,0.3,0.6,0.5]);
        clf;
        pl = plot(l_all,var_all,'.-','LineWidth',2);
        grid on;
        xlabel('$\lambda$','interpreter','latex');
        ylabel('$v_{i}$','interpreter','latex');
        xlim([max([l_lu(1),l_max(1)-dl0]),min([l_lu(2),l_max(2)+dl0])]);
        drawnow;
        pl_info = struct('fig',fig,'pl',pl);
    elseif bif_flag == -1
            %% final change in live plot
        for k = 1:nv
           row = dataTipTextRow('Step',0:l_start,'%d');
           pl_info.pl(k).DataTipTemplate.DataTipRows(end+1) = row;
        end
    else
        set(0, 'currentfigure', pl_info.fig);
        newXData = cell(nv, 1);
        newYData = cell(nv, 1);
        for k = 1:nv
            newXData{k,1} = l_all;
            newYData{k,1} = var_all(k,:);
        end
        set(pl_info.pl, {'XData'}, newXData, {'YData'},  newYData);
        if ~Opt.unique && ison(Opt.bifurcation) && bif_flag 
           if ~isempty(bif)
               hold on;
               if bif(2,end) == 0
                   plot(l_all(end),var_all(:,end),'bo','LineWidth',2);
               elseif bif(2,end) == 1
                   plot(l_all(end),var_all(:,end),'rx','LineWidth',2);
               end
               hold off;
           end
        end
        xlim([max([l_lu(1),l_max(1)-dl0]),min([l_lu(2),l_max(2)+dl0])]);
        drawnow;
    end
end

