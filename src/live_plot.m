%% path continuation - live_plot
%
%   Institute of Dynamics and Vibration Research
%   Leibniz University Hannover
%   30.09.2020 - Tido Kubatschek
%
function [pl_info,Opt] = live_plot(Opt, nv, l_start, l_end, l_all, var_all, s_all, ds, dsim1, iterations, loop_counter, fun_predictor, s_predictor, pl_info, bif_flag, bif)
    l_lu = [min([l_start,l_end]),max([l_start,l_end])];
    l_max = [min(l_all),max(l_all)];
    dl0 = abs(l_end-l_start)*0.05;
    num_pl = numel(Opt.plot_vars_index);
    if Opt.bifurcation.trace
        Opt.live_plot_fig = NaN;
    end
    if Opt.plot.basic
        if length(l_all) == 1
            
            if isnan(Opt.live_plot_fig) % test for existing figure to plot in
                %% delete used tags
                %
                if ~isempty(findobj('Tag', 'upperleft'))
                    set(findobj('Tag', 'upperleft'),'Tag','');
                end
                if ~isempty(findobj('Tag', 'upperright'))
                    set(findobj('Tag', 'upperright'),'Tag','');
                end
                if ~isempty(findobj('Tag', 'lowerleft'))
                    set(findobj('Tag', 'lowerleft'),'Tag','');
                end
                if ~isempty(findobj('Tag', 'lowerright'))
                    set(findobj('Tag', 'lowerright'),'Tag','');
                end
                %
                %% create new fig
                %
                fig = figure('units', 'normalized', 'position', [0.2,0.3,0.6,0.5]);
                clf;
            else
                fig = figure(Opt.live_plot_fig); % use existing fig
                hold on; % for new plot
            end
            
            
            %% prepare colors 
            colors = cell(num_pl,1);        
            for k = 1:num_pl
                colors(k) = {get_RGB(k,num_pl,1)};
            end
            
            %% create plot with colors
            pl = plot(l_all,var_all(Opt.plot_vars_index,:),'.-','LineWidth',2);
            set(pl, {'Color'}, colors);
            
            if isnan(Opt.live_plot_fig) || ~Opt.bifurcation.trace % test for existing figure to plot in, there must be no new labels or grid
                grid on;
                xlabel('$\lambda$','interpreter','latex');
                ylabel('$v_{i}$','interpreter','latex');
                xlim([max([l_lu(1),l_max(1)-dl0]),min([l_lu(2),l_max(2)+dl0])]);
            end
            pl_info = struct('fig',fig,'pl',pl);         
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
        else
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
                       plot(l_all(bif(1,end)),var_all(Opt.plot_vars_index,bif(1,end)),'ro','LineWidth',2);
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
        end
    elseif Opt.plot.detail
        %% find interesting plot variable
        if isnan(Opt.plot_var_of_interest)
            %% find most changing
            n_change = 20;
            [~,interesting] = max(mean(abs(diff(var_all(Opt.plot_vars_index,:),1,2))')');
        else
            %% take input by user
            interesting = Opt.plot_var_of_interest;
        end
        %
        %% get predictor and corrector
        [lco, vco] = draw_corrector(var_all, l_all, dsim1, Opt);
        l_pred = lco{1};
        l_cor = lco{2};
        v_pred = vco{1};
        v_cor = vco{2};
        %
        if length(l_all) == 1
            if isnan(Opt.live_plot_fig) % test for existing figure to plot in
                fig = figure('units', 'normalized', 'position', [0.05,0.1,0.9,0.8]); % create new fig
                clf;
            else
                fig = figure(Opt.live_plot_fig); % use existing fig
                %hold on; % for new plot
            end
            
            colors = cell(num_pl,1);        
            for k = 1:num_pl
                colors(k) = {get_RGB(k,num_pl,1)};
            end
            %% creater upper subplot
            %
            if isnan(Opt.live_plot_fig) || ~Opt.bifurcation.trace
                hs1 = subplot(2,3,1:2);
                pl = plot(l_all,var_all(Opt.plot_vars_index,:),'.-','LineWidth',2);
                set(pl, {'Color'}, colors);
                grid on;
                title('path continuation','interpreter','latex');
                xlabel('$\lambda$','interpreter','latex');
                ylabel('$v_{i}$','interpreter','latex');
                xlim([max([l_lu(1),l_max(1)-dl0]),min([l_lu(2),l_max(2)+dl0])]);
                set(hs1,'Tag','upperleft');
            else
                hs1 = findobj('Tag', 'upperleft');
                set(hs1, 'NextPlot', 'add');
                pl = plot(hs1, l_all,var_all(Opt.plot_vars_index,:),'.-','LineWidth',2);
                set(pl, {'Color'}, colors);
            end
            %
            %% create lower subplots
            %
            % first one: showing iteration steps
            %
            if isnan(Opt.live_plot_fig) || ~Opt.bifurcation.trace
                hs2 = subplot(2,3,3);
                pl_it = plot(loop_counter, iterations, 'r-x', 'LineWidth', 2);
                grid on;
                title('needed iterations per loop','interpreter','latex');
                xlabel('loop counter','interpreter','latex');
                ylabel('iterations','interpreter','latex');
                set(hs2,'Tag','upperright');
            else
                hs2 = findobj('Tag', 'upperright');
                set(hs2, 'NextPlot', 'replacechildren'); %% maybe renew ???
                pl_it = plot(hs2, loop_counter, iterations, 'r-x', 'LineWidth', 2);
            end
            %
            % second one: showing most changing variable in detail
            %
            if isnan(Opt.live_plot_fig) || ~Opt.bifurcation.trace
                hs3 = subplot(2,3,4:5);
                pl_det = plot(l_all, var_all(interesting,:),'LineWidth', 2);
                set(pl_det, {'Color'}, colors(interesting));
                grid on;
                msg = ['interesting variable: $v_',num2str(interesting),'$'];
                if isnan(Opt.plot_var_of_interest) msg = [msg, ' (most changing)']; end
                title(msg,'interpreter','latex');
                xlabel('$\lambda$','interpreter','latex');
                ylabel(['$v_',num2str(interesting),'$'],'interpreter','latex');
                xlim([l_all-2*dsim1, l_all+2*dsim1]);
                ylim([var_all(interesting)-2*dsim1, var_all(interesting)+2*dsim1]);
                set(hs3,'Tag','lowerleft');
            else
                hs3 = findobj('Tag', 'lowerleft');
                set(hs3, 'NextPlot', 'add');
                pl_det = plot(hs3, l_all, var_all(interesting,:),'LineWidth', 2);
                set(pl_det, {'Color'}, colors(interesting));
            end
            %
            % third one: showing s_all over loop_counter
            %
            if isnan(Opt.live_plot_fig) || ~Opt.bifurcation.trace
                hs4 = subplot(2,3,6);
                pl_s = plot(loop_counter, s_all,'LineWidth', 2, 'Color', 'r');
                grid on;
                title('arc length $s$','interpreter','latex');
                xlabel('loop counter','interpreter','latex');
                ylabel('$s_{\mathrm{all}}$','interpreter','latex');
                set(hs4,'Tag','lowerright');
            else
                hs4 = findobj('Tag', 'lowerright');
                set(hs4, 'NextPlot', 'add');
                pl_s = plot(hs4, loop_counter, s_all,'LineWidth', 2, 'Color', 'r');
            end
            %
            pl_info = struct('fig',fig,'pl',pl,'pl_it',pl_it,'pl_det',pl_det,'pl_s',pl_s, 'pl_pred', [], 'pl_cor', []);            
            if isnan(Opt.live_plot_fig)
                Opt.live_plot_fig = fig.Number; % reference to existing fig
            end
        elseif bif_flag == -1
            %% final change in live plot
            %
            set(0, 'currentfigure', pl_info.fig);
            %
            %% upper subplot
            %
            subplot(2,3,1:2);
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
            %% lower subplots
            %
            % first one: showing iteration steps
            %
               % no update needed
            %
            % second one: showing most changing variable in detail
            %
            subplot(2,3,4:5)
            pl_info.pl_det.XData = l_all;
            pl_info.pl_det.YData = var_all(interesting,:);
            set(pl_info.pl_det, {'Color'}, {get_RGB(interesting,num_pl,1)});
            msg = ['interesting variable: $v_',num2str(interesting),'$'];
            if isnan(Opt.plot_var_of_interest) msg = [msg, ' (most changing)']; end
            title(msg,'interpreter','latex');
            %title(['most changing variable: $v_',num2str(interesting),'$'],'interpreter','latex');
            ylabel(['$v_',num2str(interesting),'$'],'interpreter','latex');
            %
            xlim([l_all(end-1)-2*dsim1, l_all(end-1)+2*dsim1]);
            ylim([var_all(interesting,end-1)-2*dsim1, var_all(interesting,end-1)+2*dsim1]);
            hold on;
            pl_info.pl_pred = plot(l_pred(interesting,:), v_pred(interesting,:), 'r--', 'LineWidth', 1);
            pl_info.pl_cor = plot(l_cor(interesting,:), v_cor(interesting,:), 'r', 'LineWidth', 1);
            hold off;
            % third one: s_all über loop_counter
            %
                % no update needed
            %
        else
            set(0, 'currentfigure', pl_info.fig);
            %% creater upper subplot
            %
            subplot(2,3,1:2);
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
                       plot(l_all(bif(1,end)),var_all(Opt.plot_vars_index,bif(1,end)),'ro','LineWidth',2);
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
            %
            %% creat lower subplots
            %
            % first one: showing iteration steps
            %
            subplot 233
            %% save plot data
            %
            pl_info.pl_it.XData = [pl_info.pl_it.XData, loop_counter];
            pl_info.pl_it.YData = [pl_info.pl_it.YData, iterations];
            xl = xlim();
            xlim([xl(1), loop_counter]);
            %
            %
            % second one: showing most changing variable in detail
            %
            subplot(2,3,4:5);
            pl_info.pl_det.XData = l_all;
            pl_info.pl_det.YData = var_all(interesting,:);
            set(pl_info.pl_det, {'Color'}, {get_RGB(interesting,num_pl,1)});
            msg = ['interesting variable: $v_',num2str(interesting),'$'];
            if isnan(Opt.plot_var_of_interest) msg = [msg, ' (most changing)']; end
            title(msg,'interpreter','latex');
            %title(['most changing variable: $v_',num2str(interesting),'$'],'interpreter','latex');
            ylabel(['$v_',num2str(interesting),'$'],'interpreter','latex');
            %
            xlim([l_all(end-1)-2*dsim1, l_all(end-1)+2*dsim1]);
            ylim([var_all(interesting,end-1)-2*dsim1, var_all(interesting,end-1)+2*dsim1]);
            hold on;
            delete(pl_info.pl_pred);
            delete(pl_info.pl_cor);
            pl_info.pl_pred = plot(l_pred(interesting,:), v_pred(interesting,:), 'r--', 'LineWidth', 1);
            pl_info.pl_cor = plot(l_cor(interesting,:), v_cor(interesting,:), 'r', 'LineWidth', 1);
            hold off;
            %
            % third one: s_all over loop_counter
            %
            subplot 236
            pl_info.pl_s.XData = [pl_info.pl_s.XData, loop_counter];
            pl_info.pl_s.YData = [pl_info.pl_s.YData, s_all(end)];
            xlim([xl(1), loop_counter]);
            %
            %            
        end
    else
        error('No such plot method!');
    end
    drawnow limitrate;
    %
    %% Stop to show predictor and corrector
    %
    if islogical(Opt.plot_pause) && Opt.plot_pause
        fprintf('Press any key to continue or press Ctrl+c to stop...');
        pause
        fprintf(repmat('\b',1,52));
    elseif numel(l_all) >= Opt.plot_pause && ~islogical(Opt.plot_pause)
        fprintf('Press any key to continue or press Ctrl+c to stop...');
        pause
        fprintf(repmat('\b',1,52));
    end
    %
end