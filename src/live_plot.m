%% path continuation - live_plot
%
%   Institute of Dynamics and Vibration Research
%   Leibniz University Hannover
%   30.09.2020 - Tido Kubatschek
%
function [pl_info,Opt] = live_plot(Opt, nv, l_start, l_end, l_all, var_all, ds, dsim1, iterations, loop_counter, fun_predictor, s_predictor, pl_info, bif_flag, bif)
    l_lu = [min([l_start,l_end]),max([l_start,l_end])];
    l_max = [min(l_all),max(l_all)];
    dl0 = abs(l_end-l_start)*0.05;
    
    
    %% find most changing
    most_changing = 2;
    
    %%
    if length(l_all) == 1
        
        if isnan(Opt.live_plot_fig) % test for existing figure to plot in
            if Opt.plot.basic
                fig = figure('units', 'normalized', 'position', [0.2,0.3,0.6,0.5]); % create new fig
            elseif Opt.plot.detail
                fig = figure('units', 'normalized', 'position', [0.2,0.1,0.6,0.7]); % create new fig
            end
            clf;
        else
            fig = figure(Opt.live_plot_fig); % use existing fig
            hold on; % for new plot
        end
        
        num_pl = numel(var_all(:,1));
        colors = cell(num_pl,1);        
        for k = 1:num_pl
            colors(k) = {get_RGB(k,num_pl,1)};
        end
        
        if Opt.plot.basic
            pl = plot(l_all,var_all(Opt.plot_vars_index,:),'.-','LineWidth',2);
            set(pl, {'Color'}, colors);
            pl_info = struct('fig',fig,'pl',pl,'plal1',[],'plal2',[],'plpr',[]);
            if isnan(Opt.live_plot_fig) % test for existing figure to plot in, there must be no new labels or grid
                grid on;
                xlabel('$\lambda$','interpreter','latex');
                ylabel('$v_{i}$','interpreter','latex');
                xlim([max([l_lu(1),l_max(1)-dl0]),min([l_lu(2),l_max(2)+dl0])]);
            end
            
        elseif Opt.plot.detail
            %% creater upper subplot
            %
            if isnan(Opt.live_plot_fig)
                hs = subplot(2,3,1:3);
                pl = plot(l_all,var_all(Opt.plot_vars_index,:),'.-','LineWidth',2);
                set(pl, {'Color'}, colors);
                grid on;
                xlabel('$\lambda$','interpreter','latex');
                ylabel('$v_{i}$','interpreter','latex');
                xlim([max([l_lu(1),l_max(1)-dl0]),min([l_lu(2),l_max(2)+dl0])]);
                set(hs,'Tag','upper');
            else
                hs1 = findobj('Tag', 'upper');
                set(hs1, 'NextPlot', 'add');
                pl = plot(hs1, l_all,var_all(Opt.plot_vars_index,:),'.-','LineWidth',2);
                set(pl, {'Color'}, colors);
            end
            %
            %% creat lower subplots
            %
            % first one: showing iteration steps
            %
            if isnan(Opt.live_plot_fig)
                hs2 = subplot(2,3,4);
                plal1 = plot(loop_counter, iterations, 'r--.', 'LineWidth', 2);
                grid on;
                xlabel('loop counter','interpreter','latex');
                ylabel('iterations','interpreter','latex');
                set(hs2,'Tag','left');
            else
                hs2 = findobj('Tag', 'left');
                set(hs2, 'NextPlot', 'add');
                plal1 = plot(hs2, loop_counter, iterations, '--.', 'LineWidth', 2);
            end
            %
            % second one: showing most changing variable in detail
            %
            if isnan(Opt.live_plot_fig)
                hs3 = subplot(2,3,5);
                plal2 = plot(l_all, var_all(most_changing,:));
                grid on;
                xlabel('$\lambda$','interpreter','latex');
                ylabel(['$v_', num2str(2),'$'],'interpreter','latex');
                xlim([max([l_lu(1),l_max(1)-dl0]),min([l_lu(2),l_max(2)+dl0])]);
                set(hs3,'Tag','center');
            else
                hs3 = findobj('Tag', 'center');
                set(hs3, 'NextPlot', 'add');
                plal2 = plot(hs3, l_all, var_all(most_changing,:));
            end
            %
            % third one: showing predictor
            %
            if isnan(Opt.live_plot_fig)
                hs4 = subplot(2,3,6);
                plpr = plot(l_all, var_all(most_changing,:));
                grid on;
                xlabel('$\lambda$','interpreter','latex');
                ylabel('$v_{i}$','interpreter','latex');
                xlim([max([l_lu(1),l_max(1)-dl0]),min([l_lu(2),l_max(2)+dl0])]);
                set(hs4,'Tag','right');
            else
                hs4 = findobj('Tag', 'right');
                set(hs4, 'NextPlot', 'add');
                plpr = plot(hs4, l_all, var_all(most_changing,:));
            end
            %
            pl_info = struct('fig',fig,'pl',pl,'plal1',plal1,'plal2',plal2,'plpr',plpr);
        else
            error('No such plot method!');
        end
        
%         drawnow;
        
        if isnan(Opt.live_plot_fig)
            Opt.live_plot_fig = fig.Number; % reference to existing fig
        end
        
    elseif bif_flag == -1
        %% final change in live plot
        %
        set(0, 'currentfigure', pl_info.fig);
        %
        %% save plot data
        %
        newXData = cell(length(Opt.plot_vars_index), 1);
        newYData = cell(length(Opt.plot_vars_index), 1);
        for k = 1:length(Opt.plot_vars_index)
            newXData{k,1} = l_all;
            newYData{k,1} = var_all(Opt.plot_vars_index(k),:);
        end
        %
        %% plot new plot Data
        %
        set(pl_info.pl, {'XData'}, newXData, {'YData'},  newYData);
        %
        %% add third information to plot (current step)
        for k = 1:length(Opt.plot_vars_index)
            row = dataTipTextRow('Step',0:(length(l_all)-1),'%d');
            pl_info.pl(k).DataTipTemplate.DataTipRows(end+1) = row;
        end
        %
        %% adjust x axis
        xlim([max([l_lu(1),l_max(1)-dl0]),min([l_lu(2),l_max(2)+dl0])]);
        %
        drawnow;
    else
        if Opt.plot.basic
            set(0, 'currentfigure', pl_info.fig);
            %% save plot data
            %
            newXData = cell(length(Opt.plot_vars_index), 1);
            newYData = cell(length(Opt.plot_vars_index), 1);
            for k = 1:length(Opt.plot_vars_index)
                newXData{k,1} = l_all;
                newYData{k,1} = var_all(Opt.plot_vars_index(k),:);
            end
            %
            %% plot new plot Data
            %
            set(pl_info.pl, {'XData'}, newXData, {'YData'},  newYData);
            %
            %% mark bifurcation points
            %
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
            %
            %% adjust x axis
            %
            xlim([max([l_lu(1),l_max(1)-dl0]),min([l_lu(2),l_max(2)+dl0])]);
            drawnow;    
        elseif Opt.plot.detail
            set(0, 'currentfigure', pl_info.fig);
            %% creater upper subplot
            %
            subplot(2,3,1:3)
            %% save plot data
            %
            newXData = cell(length(Opt.plot_vars_index), 1);
            newYData = cell(length(Opt.plot_vars_index), 1);
            for k = 1:length(Opt.plot_vars_index)
                newXData{k,1} = l_all;
                newYData{k,1} = var_all(Opt.plot_vars_index(k),:);
            end
            %
            %% plot new plot Data
            %
            set(pl_info.pl, {'XData'}, newXData, {'YData'},  newYData);
            %
            %% mark bifurcation points
            %
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
            %
            %% adjust x axis
            %
            xlim([max([l_lu(1),l_max(1)-dl0]),min([l_lu(2),l_max(2)+dl0])]);
            drawnow;
            %
            %% creat lower subplots
            %
            % first one: showing iteration steps
            %
            subplot 234
            %% save plot data
            %
            pl_info.plal1.XData = [pl_info.plal1.XData, loop_counter];
            pl_info.plal1.YData = [pl_info.plal1.YData, iterations];
            %
            drawnow;
            %
            % second one: showing most changing variable in detail
            %
            subplot 235
            pl_info.plal2.XData = l_all;
            pl_info.plal2.YData = var_all(most_changing,:);
            %
            xlim([max([l_lu(1),l_max(1)-dl0]),min([l_lu(2),l_max(2)+dl0])]);
            drawnow;
            %
            % third one: showing predictor
            %
            subplot 236
            pl_info.plpr.XData = l_all;
            pl_info.plpr.YData = var_all(most_changing,:);
            %
            xlim([max([l_lu(1),l_max(1)-dl0]),min([l_lu(2),l_max(2)+dl0])]);
            drawnow;
            %
        else
            error('No such plot method!')
        end
    end
    
    %% detail, oben für gesamtplot trotzdem prüfen und in subplot schreiben
end

