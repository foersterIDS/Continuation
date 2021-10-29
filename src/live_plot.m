%% path continuation - live_plot
%
%   Institute of Dynamics and Vibration Research
%   Leibniz University Hannover
%   30.09.2020 - Tido Kubatschek
%
function [Plot,Opt] = live_plot(Opt, Info, Path, ds, dsim1, iterations, Counter, fun_predictor, s_predictor, Plot, Bifurcation)
    l_lu = [min([Info.l_start,Info.l_end]),max([Info.l_start,Info.l_end])];
    l_max = [min(Path.l_all),max(Path.l_all)];
    dl0 = abs(Info.l_end-Info.l_start)*0.05;
    num_pl = numel(Opt.plot_vars_index);
%     if Opt.bifurcation.trace
%         Opt.live_plot_fig = NaN;
%     end
    if (Opt.plot.basic || Opt.plot.semilogx || Opt.plot.semilogy || Opt.plot.loglog)
        if length(Path.l_all) == 1
            
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
            pl = plot(Path.l_all,Path.var_all(Opt.plot_vars_index,:),'-','LineWidth',2);
            set(pl, {'Color'}, colors); hold on;
            pl_curr = plot(Path.l_all,Path.var_all(Opt.plot_vars_index,:),'*','LineWidth',2);
            set(pl_curr, {'Color'}, colors); hold off;
            if Opt.plot.semilogx || Opt.plot.loglog
                set(gca, 'XScale', 'log')
                % check whether all lambda are positive
                if sum(Path.l_all<0)
                    error('To use semilogx-/loglog-scale lamdba must not contain negative values!');
                end
            end
            if Opt.plot.semilogy || Opt.plot.loglog
                set(gca, 'YScale', 'log')
                % check whether all variables are positive
                if sum(Path.var_all<0)
                    error('To use semilogy-/loglog-scale variables must not contain any negative values! Consinder using the abs()-function.')
                end
            end           
            
            if isnan(Opt.live_plot_fig) || ~Opt.bifurcation.trace % test for existing figure to plot in, there must be no new labels or grid
                grid on;
                xlabel('$\lambda$','interpreter','latex');
                ylabel('$v_{i}$','interpreter','latex');
                xlim([max([l_lu(1),l_max(1)-dl0]),min([l_lu(2),l_max(2)+dl0])]);
            end
            Plot.fig = fig;
            Plot.pl = pl;
            Plot.pl_curr = pl_curr;
            if isnan(Opt.live_plot_fig)
                Opt.live_plot_fig = fig.Number; % reference to existing fig
            end
            
        elseif Bifurcation.flag == -1
            %% final change in live plot
            %
            set(0, 'currentfigure', Plot.fig);
            %
            %% save plot data
            %
            newXData = cell(length(Opt.plot_vars_index), 1);
            newYData = cell(length(Opt.plot_vars_index), 1);
            newXData_curr = cell(length(Opt.plot_vars_index), 1);
            newYData_curr = cell(length(Opt.plot_vars_index), 1);
            for k = 1:length(Opt.plot_vars_index)
                newXData{k,1} = Path.l_all;
                newYData{k,1} = Path.var_all(Opt.plot_vars_index(k),:);
%                 newXData_curr{k,1} = Path.l_all(end);
%                 newYData_curr{k,1} = Path.var_all(Opt.plot_vars_index(k),end);
            end
            %
            %% plot new plot Data
            %
            set(Plot.pl, {'XData'}, newXData, {'YData'},  newYData);
%             set(Plot.pl_curr, {'XData'}, newXData_curr, {'YData'},  newYData_curr);
            %
            %% add third information to plot (current step)
            for k = 1:length(Opt.plot_vars_index)
                row = dataTipTextRow('Step',0:(length(Path.l_all)-1),'%d');
                Plot.pl(k).DataTipTemplate.DataTipRows(end+1) = row;
            end
            %
            %% adjust x axis
            if Opt.bifurcation.trace
                axis([Info.l_start, Info.l_end, min(min(Path.var_all)), max(max(Path.var_all))]);
            else
                xlim([max([l_lu(1),l_max(1)-dl0]),min([l_lu(2),l_max(2)+dl0])]);
            end
        else
            set(0, 'currentfigure', Plot.fig);
            %% save plot data
            %
            newXData = cell(length(Opt.plot_vars_index), 1);
            newYData = cell(length(Opt.plot_vars_index), 1);
            newXData_curr = cell(length(Opt.plot_vars_index), 1);
            newYData_curr = cell(length(Opt.plot_vars_index), 1);
            for k = 1:length(Opt.plot_vars_index)
                newXData{k,1} = Path.l_all;
                newYData{k,1} = Path.var_all(Opt.plot_vars_index(k),:);
                newXData_curr{k,1} = Path.l_all(end);
                newYData_curr{k,1} = Path.var_all(Opt.plot_vars_index(k),end);
            end
            %
            %% plot new plot Data
            %
            set(Plot.pl, {'XData'}, newXData, {'YData'},  newYData);
            set(Plot.pl_curr, {'XData'}, newXData_curr, {'YData'},  newYData_curr);
            %
            %% mark bifurcation points
            %
            if ison(Opt.bifurcation) && Bifurcation.flag 
               if ~isempty(Bifurcation.bif)
                   hold on;
                   if Bifurcation.bif(2,end) == 0 % brach point Bifurcation.bif
                       plot(Path.l_all(Bifurcation.bif(1,end)),Path.var_all(Opt.plot_vars_index,Bifurcation.bif(1,end)),'rx','LineWidth',2);
                   elseif Bifurcation.bif(2,end) == 1 % fold Bifurcation.bif
                       plot(Path.l_all(Bifurcation.bif(1,end)),Path.var_all(Opt.plot_vars_index,Bifurcation.bif(1,end)),'ro','LineWidth',2);
                   elseif isnan(Bifurcation.bif(2,end)) % Bifurcation.bif, but no further information
                       plot(Path.l_all(Bifurcation.bif(1,end)),Path.var_all(Opt.plot_vars_index,Bifurcation.bif(1,end)),'rs','LineWidth',2);
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
            [~,interesting] = max(mean(abs(diff(Path.var_all(Opt.plot_vars_index,:),1,2))')');
        else
            %% take input by user
            interesting = Opt.plot_var_of_interest;
        end
        %
        %% get corrector
        [lco, vco] = draw_corrector(Path, dsim1, Opt);
        l_cor_assist = lco{1};
        l_cor = lco{2};
        v_cor_assist = vco{1};
        v_cor = vco{2};
        %
        %% get predictor
        if nargin>11
            s_pre = linspace(0,s_predictor,50);
            x_pre = fun_predictor(s_pre);
            l_pre = kron(ones(Info.nv,1),x_pre(end,:));
            v_pre = x_pre(1:(end-1),:);
        else
            l_pre = NaN(Info.nv,1);
            v_pre = NaN(Info.nv,1);
        end
        %
        if length(Path.l_all) == 1
            if isnan(Opt.live_plot_fig) % test for existing figure to plot in
                fig = figure('units', 'normalized', 'position', [0.05,0.1,0.9,0.8]); % create new fig
                clf;
            else
                fig = figure(Opt.live_plot_fig); % use existing fig
                hold on; % for new plot
            end
            
            colors = cell(num_pl,1);        
            for k = 1:num_pl
                colors(k) = {get_RGB(k,num_pl,1)};
            end
            %% creater upper subplot
            %
            if isnan(Opt.live_plot_fig) || ~Opt.bifurcation.trace
                hs1 = subplot(2,3,1:2);
                hold on;
                pl = plot(Path.l_all,Path.var_all(Opt.plot_vars_index,:),'-','LineWidth',2);
                set(pl, {'Color'}, colors);
                pl_curr = plot(Path.l_all,Path.var_all(Opt.plot_vars_index,:),'*','LineWidth',2);
                set(pl_curr, {'Color'}, colors);
                grid on;
                title('path continuation','interpreter','latex');
                xlabel('$\lambda$','interpreter','latex');
                ylabel('$v_{i}$','interpreter','latex');
                xlim([max([l_lu(1),l_max(1)-dl0]),min([l_lu(2),l_max(2)+dl0])]);
                set(hs1,'Tag','upperleft');
            else
                hs1 = findobj('Tag', 'upperleft');
                set(hs1, 'NextPlot', 'add');
                pl = plot(hs1, Path.l_all,Path.var_all(Opt.plot_vars_index,:),'-','LineWidth',2);
                set(pl, {'Color'}, colors);
                pl_curr = plot(Path.l_all,Path.var_all(Opt.plot_vars_index,:),'*','LineWidth',2);
                set(pl_curr, {'Color'}, colors);
            end
            %
            % second one: showing iteration steps
            %
            if isnan(Opt.live_plot_fig) || ~Opt.bifurcation.trace
                hs2 = subplot(2,3,3);
                pl_it = plot(Counter.loop, iterations, 'r-x', 'LineWidth', 2);
                grid on;
                title('needed iterations per loop','interpreter','latex');
                xlabel('loop counter','interpreter','latex');
                ylabel('iterations','interpreter','latex');
                set(hs2,'Tag','upperright');
            else
                hs2 = findobj('Tag', 'upperright');
                set(hs2, 'NextPlot', 'replacechildren'); %% maybe renew ???
                pl_it = plot(hs2, Counter.loop, iterations, 'r-x', 'LineWidth', 2);
            end
            %
            %% create lower subplots
            %
            % first one: showing most changing variable in detail
            %
            if isnan(Opt.live_plot_fig) || ~Opt.bifurcation.trace
                hs3 = subplot(2,3,4:5);
                pl_det = plot(Path.l_all, Path.var_all(interesting,:),'LineWidth', 2);
                set(pl_det, {'Color'}, colors(interesting));
                grid on;
                msg = ['interesting variable: $v_{',num2str(interesting),'}$'];
                if isnan(Opt.plot_var_of_interest) msg = [msg, ' (most changing)']; end
                title(msg,'interpreter','latex');
                xlabel('$\lambda$','interpreter','latex');
                ylabel(['$v_{',num2str(interesting),'}$'],'interpreter','latex');
                xlim([Path.l_all-2*dsim1, Path.l_all+2*dsim1]);
                v_diff = abs(v_pre(interesting,end) - v_pre(interesting,1));
                if 1.1*v_diff > 2*dsim1
                    ylim([Path.var_all(interesting)-v_diff, Path.var_all(interesting)+v_diff]);
                else
                    ylim([Path.var_all(interesting)-2*dsim1, Path.var_all(interesting)+2*dsim1]);
                end
                hold on;
                pl_cor_assist = plot(l_cor_assist(interesting,:), v_cor_assist(interesting,:), 'r--', 'LineWidth', 1);
                pl_cor = plot(l_cor(interesting,:), v_cor(interesting,:), 'r', 'LineWidth', 1);
                pl_pre = plot(l_pre(interesting,:), v_pre(interesting,:), 'k-.', 'LineWidth', 1);
%                 hold off;
                set(hs3,'Tag','lowerleft');
            else
                hs3 = findobj('Tag', 'lowerleft');
                set(hs3, 'NextPlot', 'add');
                pl_det = plot(hs3, Path.l_all, Path.var_all(interesting,:),'LineWidth', 2);
                set(pl_det, {'Color'}, colors(interesting));
            end
            %
            % second one: showing Path.s_all over Counter.loop
            %
            if isnan(Opt.live_plot_fig) || ~Opt.bifurcation.trace
                hs4 = subplot(2,3,6);
                hold off;
                pl_s = plot(Counter.loop, Path.s_all,'LineWidth', 2, 'Color', 'r');
                grid on;
                title('arc length $s$','interpreter','latex');
                xlabel('loop counter','interpreter','latex');
                ylabel('$s_{\mathrm{all}}$','interpreter','latex');
                set(hs4,'Tag','lowerright');
            else
                hs4 = findobj('Tag', 'lowerright');
                set(hs4, 'NextPlot', 'add');
                pl_s = plot(hs4, Counter.loop, Path.s_all,'LineWidth', 2, 'Color', 'r');
            end
            %
            Plot.fig = fig; Plot.pl = pl; Plot.pl_it = pl_it;
            Plot.pl_det = pl_det; Plot.pl_s = pl_s; Plot.pl_cor_assist = pl_cor_assist;
            Plot.pl_cor = pl_cor; Plot.pl_pre = pl_pre;
            Plot.pl_curr = pl_curr;
            if isnan(Opt.live_plot_fig)
                Opt.live_plot_fig = fig.Number; % reference to existing fig
            end
        elseif Bifurcation.flag == -1
            %% final change in live plot
            %
            set(0, 'currentfigure', Plot.fig);
            %
            %% calc limits
            dl = abs(Path.l_all(end)-Path.l_all(end-1));
            dv = abs(Path.var_all(interesting,end)-Path.var_all(interesting,end-1));
            %
            %% upper subplots
            %
            % first one
            %
            subplot(2,3,1:2);
            %
            %% save plot data
            %
            newXData = cell(length(Opt.plot_vars_index), 1);
            newYData = cell(length(Opt.plot_vars_index), 1);
%             newXData_curr = cell(length(Opt.plot_vars_index), 1);
%             newYData_curr = cell(length(Opt.plot_vars_index), 1);
            for k = 1:length(Opt.plot_vars_index)
                newXData{k,1} = Path.l_all;
                newYData{k,1} = Path.var_all(Opt.plot_vars_index(k),:);
%                 newXData_curr{k,1} = Path.l_all(end);
%                 newYData_curr{k,1} = Path.var_all(Opt.plot_vars_index(k),end);
            end
            %
            %% plot new plot Data
            %
            set(Plot.pl, {'XData'}, newXData, {'YData'},  newYData);
%             set(Plot.pl_curr, {'XData'}, newXData_curr, {'YData'},  newYData_curr);
            %
            %% add third information to plot (current step)
            for k = 1:length(Opt.plot_vars_index)
                row = dataTipTextRow('Step',0:(length(Path.l_all)-1),'%d');
                Plot.pl(k).DataTipTemplate.DataTipRows(end+1) = row;
            end
            %
            %% adjust x axis
            if Opt.bifurcation.trace
                axis([Info.l_start, Info.l_end, min(min(Path.var_all)), max(max(Path.var_all))]);
            else
                xlim([max([l_lu(1),l_max(1)-dl0]),min([l_lu(2),l_max(2)+dl0])]);
            end
            %
            % second one: showing iteration steps
            %
               % no update needed
            %
            %% lower subplots
            % first one: showing most changing variable in detail
            %
            subplot(2,3,4:5)
            Plot.pl_det.XData = Path.l_all;
            Plot.pl_det.YData = Path.var_all(interesting,:);
            set(Plot.pl_det, {'Color'}, {get_RGB(interesting,num_pl,1)});
            msg = ['interesting variable: $v_{',num2str(interesting),'}$'];
            if isnan(Opt.plot_var_of_interest) msg = [msg, ' (most changing)']; end
            title(msg,'interpreter','latex');
            ylabel(['$v_{',num2str(interesting),'}$'],'interpreter','latex');
            %
            xlim([Path.l_all(end-1)-2*dl, Path.l_all(end-1)+2*dl]);
            v_diff = abs(v_pre(interesting,end) - v_pre(interesting,1));
            if 1.1*v_diff > 2*dv
                ylim([Path.var_all(interesting,end-1)-v_diff, Path.var_all(interesting,end-1)+v_diff]);
            else
                ylim([Path.var_all(interesting,end-1)-2*dv, Path.var_all(interesting,end-1)+2*dv]);
            end
            Plot.pl_cor_assist.XData = l_cor_assist(interesting,:);
            Plot.pl_cor_assist.YData = v_cor_assist(interesting,:);
            Plot.pl_cor.XData = l_cor(interesting,:);
            Plot.pl_cor.YData = v_cor(interesting,:);
            Plot.pl_pre.XData = l_pre(interesting,:);
            Plot.pl_pre.YData = v_pre(interesting,:);
            %
            % second one: Path.s_all over Counter.loop
            %
                % no update needed
            %
        else
            set(0, 'currentfigure', Plot.fig);
            %% calc limits
            dl = abs(Path.l_all(end)-Path.l_all(end-1));
            dv = abs(Path.var_all(interesting,end)-Path.var_all(interesting,end-1));
            %
            %% creater upper subplots
            %
            % first one
            %
            subplot(2,3,1:2);
            %% save plot data
            %
            newXData = cell(length(Opt.plot_vars_index), 1);
            newYData = cell(length(Opt.plot_vars_index), 1);
            newXData_curr = cell(length(Opt.plot_vars_index), 1);
            newYData_curr = cell(length(Opt.plot_vars_index), 1);
            for k = 1:length(Opt.plot_vars_index)
                newXData{k,1} = Path.l_all;
                newYData{k,1} = Path.var_all(Opt.plot_vars_index(k),:);
                newXData_curr{k,1} = Path.l_all(end);
                newYData_curr{k,1} = Path.var_all(Opt.plot_vars_index(k),end);
            end
            %
            %% plot new plot Data
            %
            set(Plot.pl, {'XData'}, newXData, {'YData'},  newYData);
            set(Plot.pl_curr, {'XData'}, newXData_curr, {'YData'},  newYData_curr);
            %
            %% mark bifurcation points
            %
            if ison(Opt.bifurcation) && Bifurcation.flag 
               if ~isempty(Bifurcation.bif)
                   hold on;
                   if Bifurcation.bif(2,end) == 0 % brach point Bifurcation.bif
                       plot(Path.l_all(Bifurcation.bif(1,end)),Path.var_all(Opt.plot_vars_index,Bifurcation.bif(1,end)),'rx','LineWidth',2);
                   elseif Bifurcation.bif(2,end) == 1 % fold Bifurcation.bif
                       plot(Path.l_all(Bifurcation.bif(1,end)),Path.var_all(Opt.plot_vars_index,Bifurcation.bif(1,end)),'ro','LineWidth',2);
                   elseif isnan(Bifurcation.bif(2,end)) % Bifurcation.bif, but no further information
                       plot(Path.l_all(Bifurcation.bif(1,end)),Path.var_all(Opt.plot_vars_index,Bifurcation.bif(1,end)),'rs','LineWidth',2);
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
            Plot.pl_it.XData = [Plot.pl_it.XData, Counter.loop];
            Plot.pl_it.YData = [Plot.pl_it.YData, iterations];
            xl = xlim();
            xlim([xl(1), Counter.loop]);
            %
            %
            % second one: showing most changing variable in detail
            %
            subplot(2,3,4:5);
            Plot.pl_det.XData = Path.l_all;
            Plot.pl_det.YData = Path.var_all(interesting,:);
            set(Plot.pl_det, {'Color'}, {get_RGB(interesting,num_pl,1)});
            msg = ['interesting variable: $v_{',num2str(interesting),'}$'];
            if isnan(Opt.plot_var_of_interest) msg = [msg, ' (most changing)']; end
            title(msg,'interpreter','latex');
            ylabel(['$v_{',num2str(interesting),'}$'],'interpreter','latex');
            %
            xlim([Path.l_all(end-1)-2*dl, Path.l_all(end-1)+2*dl]);
            v_diff = abs(v_pre(interesting,end) - v_pre(interesting,1));
            if 1.1*v_diff > 2*dv
                ylim([Path.var_all(interesting,end-1)-v_diff, Path.var_all(interesting,end-1)+v_diff]);
            else
                ylim([Path.var_all(interesting,end-1)-2*dv, Path.var_all(interesting,end-1)+2*dv]);
            end
            Plot.pl_cor_assist.XData = l_cor_assist(interesting,:);
            Plot.pl_cor_assist.YData = v_cor_assist(interesting,:);
            Plot.pl_cor.XData = l_cor(interesting,:);
            Plot.pl_cor.YData = v_cor(interesting,:);
            Plot.pl_pre.XData = l_pre(interesting,:);
            Plot.pl_pre.YData = v_pre(interesting,:);
            %
            %
            % third one: Path.s_all over Counter.loop
            %
            subplot 236
            Plot.pl_s.XData = [Plot.pl_s.XData, Counter.loop];
            Plot.pl_s.YData = [Plot.pl_s.YData, Path.s_all(end)];
            xlim([xl(1), Counter.loop]);
            %
            %            
        end
    elseif Opt.plot.three_dim
        if num_pl+1 > 3
            warning('3D plot only works with two variables.\nFirst two are selected! Consider defining plot_vars_index.');
            Opt.plot_vars_index = [1,2];
            num_pl = 1;
        end
        if length(Path.l_all) == 1
            
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
            color = get_RGB(1,num_pl,1);

            
            %% create plot with colors
            pl = plot3(Path.l_all,Path.var_all(Opt.plot_vars_index(1),:),Path.var_all(Opt.plot_vars_index(2),:),'-','LineWidth',2);
            set(pl, 'Color', color); hold on;
            pl_curr = plot3(Path.l_all,Path.var_all(Opt.plot_vars_index(1),:),Path.var_all(Opt.plot_vars_index(2),:),'*','LineWidth',2);
            set(pl_curr, 'Color', color); hold off;
            [caz,cel] = view();
            view(caz+180,cel);
            if isnan(Opt.live_plot_fig) || ~Opt.bifurcation.trace % test for existing figure to plot in, there must be no new labels or grid
                grid on;
                xlabel('$\lambda$','interpreter','latex');
                ylabel(['$v_{',num2str(Opt.plot_vars_index(1)),'}$'],'interpreter','latex');
                zlabel(['$v_{',num2str(Opt.plot_vars_index(2)),'}$'],'interpreter','latex');
                xlim([max([l_lu(1),l_max(1)-dl0]),min([l_lu(2),l_max(2)+dl0])]);
            end
            Plot.fig = fig;
            Plot.pl = pl;
            Plot.pl_curr = pl_curr;
            if isnan(Opt.live_plot_fig)
                Opt.live_plot_fig = fig.Number; % reference to existing fig
            end
            
        elseif Bifurcation.flag == -1
            %% final change in live plot
            %
            set(0, 'currentfigure', Plot.fig);
            %
            %% save plot data
            %
            newXData = Path.l_all;
            newYData = Path.var_all(Opt.plot_vars_index(1),:);
            newZData = Path.var_all(Opt.plot_vars_index(2),:);
%             newXData_curr = Path.l_all(end);
%             newYData_curr = Path.var_all(Opt.plot_vars_index(1),end);
%             newZData_curr = Path.var_all(Opt.plot_vars_index(2),end);
            %
            %% plot new plot Data
            %
            set(Plot.pl, 'XData', newXData, 'YData',  newYData, 'ZData', newZData);
%             set(Plot.pl_curr, {'XData'}, newXData_curr, {'YData'},  newYData_curr);
            %
            %% add third information to plot (current step)
            row = dataTipTextRow('Step',0:(length(Path.l_all)-1),'%d');
            Plot.pl.DataTipTemplate.DataTipRows(end+1) = row;
            %
            %% adjust x axis
            if Opt.bifurcation.trace
                axis([Info.l_start, Info.l_end, min(min(Path.var_all)), max(max(Path.var_all))]);
            else
                xlim([max([l_lu(1),l_max(1)-dl0]),min([l_lu(2),l_max(2)+dl0])]);
            end
        else
            set(0, 'currentfigure', Plot.fig);
            %% save plot data
            %
            newXData = Path.l_all;
            newYData = Path.var_all(Opt.plot_vars_index(1),:);
            newZData = Path.var_all(Opt.plot_vars_index(2),:);
            newXData_curr = Path.l_all(end);
            newYData_curr = Path.var_all(Opt.plot_vars_index(1),end);
            newZData_curr = Path.var_all(Opt.plot_vars_index(2),end);
            %
            %% plot new plot Data
            %
            set(Plot.pl, 'XData', newXData, 'YData',  newYData, 'ZData', newZData);
            set(Plot.pl_curr, 'XData', newXData_curr, 'YData',  newYData_curr, 'ZData', newZData_curr);
            %
            %% mark bifurcation points
            %
            if ison(Opt.bifurcation) && Bifurcation.flag 
               if ~isempty(Bifurcation.bif)
                   hold on;
                   if Bifurcation.bif(2,end) == 0 % brach point Bifurcation.bif
                       plot(Path.l_all(Bifurcation.bif(1,end)),Path.var_all(Opt.plot_vars_index,Bifurcation.bif(1,end)),'rx','LineWidth',2);
                   elseif Bifurcation.bif(2,end) == 1 % fold Bifurcation.bif
                       plot(Path.l_all(Bifurcation.bif(1,end)),Path.var_all(Opt.plot_vars_index,Bifurcation.bif(1,end)),'ro','LineWidth',2);
                   elseif isnan(Bifurcation.bif(2,end)) % Bifurcation.bif, but no further information
                       plot(Path.l_all(Bifurcation.bif(1,end)),Path.var_all(Opt.plot_vars_index,Bifurcation.bif(1,end)),'rs','LineWidth',2);
                   end
                   hold off;
               end
            end
            %
            %% adjust x axis
            %
            xlim([max([l_lu(1),l_max(1)-dl0]),min([l_lu(2),l_max(2)+dl0])]);        
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
    elseif numel(Path.l_all) >= Opt.plot_pause && ~islogical(Opt.plot_pause)
        fprintf('Press any key to continue or press Ctrl+c to stop...');
        pause
        fprintf(repmat('\b',1,52));
    end
    %
end