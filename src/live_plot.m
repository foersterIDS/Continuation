%% path continuation - live_plot
%
%   Institute of Dynamics and Vibration Research
%   Leibniz University Hannover
%   30.09.2020 - Tido Kubatschek
%
function [pl_info,Opt] = live_plot(Opt, nv, l_start, l_end, l_all, var_all, ds, dsim1, pl_info, bif_flag, bif)
    l_lu = [min([l_start,l_end]),max([l_start,l_end])];
    l_max = [min(l_all),max(l_all)];
    dl0 = abs(l_end-l_start)*0.05;
    
    if length(l_all) == 1
        if isnan(Opt.live_plot_fig) % test for existing figure to plot in
            fig = figure('units', 'normalized', 'position', [0.2,0.3,0.6,0.5]); % create new fig
            clf;
        else
            fig = figure(Opt.live_plot_fig); % use existing fig
            hold on; % for new plot
        end
        
        if Opt.bifurcation.mark || ~ison(Opt.bifurcation)
            pl = plot(l_all,var_all(Opt.plot_vars_index,:),'.-','LineWidth',2);
        elseif Opt.bifurcation.trace || Opt.bifurcation.determine
            pl = plot(l_all,var_all(Opt.plot_vars_index,:),'.-','LineWidth',2, 'Color', [0 0.4470 0.7410]);
        end
        
        if isnan(Opt.live_plot_fig) % test for existing figure to plot in, there must to no new labels or grid
            grid on;
            xlabel('$\lambda$','interpreter','latex');
            ylabel('$v_{i}$','interpreter','latex');
            xlim([max([l_lu(1),l_max(1)-dl0]),min([l_lu(2),l_max(2)+dl0])]);
        end
        
        drawnow;
        pl_info = struct('fig',fig,'pl',pl);
        
        if isnan(Opt.live_plot_fig)
            Opt.live_plot_fig = fig.Number; % reference to existing fig
        end
        
    elseif bif_flag == -1
        %% final change in live plot
        set(0, 'currentfigure', pl_info.fig);
        newXData = cell(length(Opt.plot_vars_index), 1);
        newYData = cell(length(Opt.plot_vars_index), 1);
        for k = 1:length(Opt.plot_vars_index)
            newXData{k,1} = l_all;
            newYData{k,1} = var_all(Opt.plot_vars_index(k),:);
        end
        set(pl_info.pl, {'XData'}, newXData, {'YData'},  newYData);
        for k = 1:length(Opt.plot_vars_index)
            row = dataTipTextRow('Step',0:(length(l_all)-1),'%d');
            pl_info.pl(k).DataTipTemplate.DataTipRows(end+1) = row;
        end
        xlim([max([l_lu(1),l_max(1)-dl0]),min([l_lu(2),l_max(2)+dl0])]);
        drawnow;
    else
        set(0, 'currentfigure', pl_info.fig);
        newXData = cell(length(Opt.plot_vars_index), 1);
        newYData = cell(length(Opt.plot_vars_index), 1);
        for k = 1:length(Opt.plot_vars_index)
            newXData{k,1} = l_all;
            newYData{k,1} = var_all(Opt.plot_vars_index(k),:);
        end
        set(pl_info.pl, {'XData'}, newXData, {'YData'},  newYData);
        if ison(Opt.bifurcation) && bif_flag 
           if ~isempty(bif)
               hold on;
               if bif(2,end) == 0
                   plot(l_all(bif(1,end)),var_all(Opt.plot_vars_index,bif(1,end)),'bo','LineWidth',2);
               elseif bif(2,end) == 1
                   plot(l_all(bif(1,end)),var_all(Opt.plot_vars_index,bif(1,end)),'rx','LineWidth',2);
               end
               hold off;
           end
        end
        xlim([max([l_lu(1),l_max(1)-dl0]),min([l_lu(2),l_max(2)+dl0])]);
        drawnow;
    end
end

