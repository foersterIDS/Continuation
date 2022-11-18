%% path continuation - plot.livePlot
%
%   Institute of Dynamics and Vibration Research
%   Leibniz University Hannover
%   30.09.2020 - Tido Kubatschek
%   07.08.2022 - Anna Lefken
%
function [Plot,Opt] = livePlot(Opt, Info, Path, ds, dsim1, iterations, Counter, funPredictor, sPredictor, Plot, Bifurcation, dpaPoints)
    lLu = [min([Info.lStart,Info.lEnd]),max([Info.lStart,Info.lEnd])];
    lMax = [min(Path.lAll),max(Path.lAll)];
    if numel(Path.lAll)>=2
        dl0 = abs(max(Path.lAll)-min(Path.lAll))*0.2;
    else
        dl0 = abs(Info.lEnd-Info.lStart);
    end
    numPl = numel(Opt.plotVarsIndex);
    if nargin<12
        dpaPoints = [];
    end
%     if Opt.bifurcation.trace
%         Opt.livePlotFig = NaN;
%     end
    if (Opt.plot.basic || Opt.plot.semilogx || Opt.plot.semilogy || Opt.plot.loglog)
        if length(Path.lAll) == 1
            
            if isnan(Opt.livePlotFig) % test for existing figure to plot in
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
                fig = figure(Opt.livePlotFig); % use existing fig
                hold on; % for new plot
            end
            %% prepare colors
            colors = cell(numPl,1);        
            for k = 1:numPl
                colors(k) = {plot.getRGB(k,numPl,1)};
            end
            %% create plot with colors
            pl = plot(Path.lAll,Path.varAll(Opt.plotVarsIndex,:),'-','LineWidth',2);
            set(pl, {'Color'}, colors); hold on;
            plCurr = plot(Path.lAll,Path.varAll(Opt.plotVarsIndex,:),'*','LineWidth',2);
            set(plCurr, {'Color'}, colors); hold off;
            if Opt.plot.semilogx || Opt.plot.loglog
                set(gca, 'XScale', 'log')
                % check whether all lambda are positive
                if sum(Path.lAll<0)
                    error('To use semilogx-/loglog-scale lamdba must not contain negative values!');
                end
            end
            if Opt.plot.semilogy || Opt.plot.loglog
                set(gca, 'YScale', 'log')
                % check whether all variables are positive
                if sum(Path.varAll<0)
                    error('To use semilogy-/loglog-scale variables must not contain any negative values! Consinder using the abs()-function.')
                end
            end           
            
            if isnan(Opt.livePlotFig) || ~Opt.bifurcation.trace % test for existing figure to plot in, there must be no new labels or grid
                grid on;
                xlabel('$\lambda$','interpreter','latex');
                ylabel('$v_{i}$','interpreter','latex');
                xlim([max([lLu(1),lMax(1)-dl0]),min([lLu(2),lMax(2)+dl0])]);
            end
            Plot.fig = fig;
            Plot.pl = pl;
            Plot.plCurr = plCurr;
            if isnan(Opt.livePlotFig)
                Opt.livePlotFig = fig.Number; % reference to existing fig
            end
            
        elseif Bifurcation.flag == -1
            %% final change in live plot
            %
            set(0, 'currentfigure', Plot.fig);
            %
            %% save plot data
            %
            newXData = cell(length(Opt.plotVarsIndex), 1);
            newYData = cell(length(Opt.plotVarsIndex), 1);
            newXDataCurr = cell(length(Opt.plotVarsIndex), 1);
            newYDataCurr = cell(length(Opt.plotVarsIndex), 1);
            for k = 1:length(Opt.plotVarsIndex)
                newXData{k,1} = Path.lAll;
                newYData{k,1} = Path.varAll(Opt.plotVarsIndex(k),:);
%                 newXDataCurr{k,1} = Path.lAll(end);
%                 newYDataCurr{k,1} = Path.varAll(Opt.plotVarsIndex(k),end);
            end
            %
            %% plot new plot Data
            %
            set(Plot.pl, {'XData'}, newXData, {'YData'},  newYData);
%             set(Plot.plCurr, {'XData'}, newXDataCurr, {'YData'},  newYDataCurr);
            %
            %% add third information to plot (current step)
            for k = 1:length(Opt.plotVarsIndex)
                row = dataTipTextRow('Step',0:(length(Path.lAll)-1),'%d');
                Plot.pl(k).DataTipTemplate.DataTipRows(end+1) = row;
            end
            %
            %% adjust x axis
            if Opt.bifurcation.trace
                axis([Info.lStart, Info.lEnd, min(min(Path.varAll)), max(max(Path.varAll))]);
            else
                xlim([max([lLu(1),lMax(1)-dl0]),min([lLu(2),lMax(2)+dl0])]);
            end
        else
            set(0, 'currentfigure', Plot.fig);
            %% save plot data
            %
            newXData = cell(length(Opt.plotVarsIndex), 1);
            newYData = cell(length(Opt.plotVarsIndex), 1);
            newXDataCurr = cell(length(Opt.plotVarsIndex), 1);
            newYDataCurr = cell(length(Opt.plotVarsIndex), 1);
            for k = 1:length(Opt.plotVarsIndex)
                newXData{k,1} = Path.lAll;
                newYData{k,1} = Path.varAll(Opt.plotVarsIndex(k),:);
                newXDataCurr{k,1} = Path.lAll(end);
                newYDataCurr{k,1} = Path.varAll(Opt.plotVarsIndex(k),end);
            end
            %
            %% plot new plot Data
            %
            set(Plot.pl, {'XData'}, newXData, {'YData'},  newYData);
            set(Plot.plCurr, {'XData'}, newXDataCurr, {'YData'},  newYDataCurr);
            %
            %% mark bifurcation points
            %
            if aux.ison(Opt.bifurcation) && Bifurcation.flag 
               if ~isempty(Bifurcation.bif)
                   hold on;
                   if Bifurcation.bif(2,end) == 0 % brach point Bifurcation.bif
                       plot(Path.lAll(Bifurcation.bif(1,end)),Path.varAll(Opt.plotVarsIndex,Bifurcation.bif(1,end)),'rx','LineWidth',2);
                   elseif Bifurcation.bif(2,end) == 1 % fold Bifurcation.bif
                       plot(Path.lAll(Bifurcation.bif(1,end)),Path.varAll(Opt.plotVarsIndex,Bifurcation.bif(1,end)),'ro','LineWidth',2);
                   elseif Bifurcation.bif(2,end) == 2 % Bifurkation from additional testfunction
                       plot(Path.lAll(Bifurcation.bif(1,end)),Path.varAll(Opt.plotVarsIndex,Bifurcation.bif(1,end)),'rd','LineWidth',2);
                   elseif isnan(Bifurcation.bif(2,end)) % Bifurcation.bif, but no further information
                       plot(Path.lAll(Bifurcation.bif(1,end)),Path.varAll(Opt.plotVarsIndex,Bifurcation.bif(1,end)),'rs','LineWidth',2);
                   end
                   hold off;
               end
            end
            %
            %% adjust x axis
            %
            xlim([max([lLu(1),lMax(1)-dl0]),min([lLu(2),lMax(2)+dl0])]);        
        end
    elseif Opt.plot.detail
        %% find interesting plot variable
        if isnan(Opt.plotVarOfInterest)
            %% find most changing
            nChange = 20;
            [~,interesting] = max(mean(abs(diff(Path.varAll(Opt.plotVarsIndex,:),1,2))')');
        else
            %% take input by user
            interesting = Opt.plotVarOfInterest;
        end
        %
        %% get corrector
        [lco, vco] = plot.drawCorrector(Path, dsim1, Opt);
        lCorAssist = lco{1};
        lCor = lco{2};
        vCorAssist = vco{1};
        vCor = vco{2};
        %
        %% get predictor
        if nargin>11
            sPre = linspace(0,sPredictor,50);
            xPre = funPredictor(sPre);
            lPre = kron(ones(Info.nv,1),xPre(end,:));
            vPre = xPre(1:(end-1),:);
        else
            lPre = NaN(Info.nv,1);
            vPre = NaN(Info.nv,1);
        end
        %
        if length(Path.lAll) == 1
            if isnan(Opt.livePlotFig) % test for existing figure to plot in
                fig = figure('units', 'normalized', 'position', [0.05,0.1,0.9,0.8]); % create new fig
                clf;
            else
                fig = figure(Opt.livePlotFig); % use existing fig
                hold on; % for new plot
            end
            
            colors = cell(numPl,1);        
            for k = 1:numPl
                colors(k) = {plot.getRGB(k,numPl,1)};
            end
            %% creater upper subplot
            %
            if isnan(Opt.livePlotFig) || ~Opt.bifurcation.trace
                hs1 = subplot(2,3,1:2);
                hold on;
                pl = plot(Path.lAll,Path.varAll(Opt.plotVarsIndex,:),'-','LineWidth',2);
                set(pl, {'Color'}, colors);
                plCurr = plot(Path.lAll,Path.varAll(Opt.plotVarsIndex,:),'*','LineWidth',2);
                set(plCurr, {'Color'}, colors);
                grid on;
                title('path continuation','interpreter','latex');
                xlabel('$\lambda$','interpreter','latex');
                ylabel('$v_{i}$','interpreter','latex');
                xlim([max([lLu(1),lMax(1)-dl0]),min([lLu(2),lMax(2)+dl0])]);
                set(hs1,'Tag','upperleft');
            else
                hs1 = findobj('Tag', 'upperleft');
                set(hs1, 'NextPlot', 'add');
                pl = plot(hs1, Path.lAll,Path.varAll(Opt.plotVarsIndex,:),'-','LineWidth',2);
                set(pl, {'Color'}, colors);
                plCurr = plot(Path.lAll,Path.varAll(Opt.plotVarsIndex,:),'*','LineWidth',2);
                set(plCurr, {'Color'}, colors);
            end
            %
            % second one: showing iteration steps
            %
            if isnan(Opt.livePlotFig) || ~Opt.bifurcation.trace
                hs2 = subplot(2,3,3);
                plIt = plot(Counter.loop, iterations, 'r-x', 'LineWidth', 2);
                grid on;
                title('needed iterations per loop','interpreter','latex');
                xlabel('loop counter','interpreter','latex');
                ylabel('iterations','interpreter','latex');
                set(hs2,'Tag','upperright');
            else
                hs2 = findobj('Tag', 'upperright');
                set(hs2, 'NextPlot', 'replacechildren'); %% maybe renew ???
                plIt = plot(hs2, Counter.loop, iterations, 'r-x', 'LineWidth', 2);
            end
            %
            %% create lower subplots
            %
            % first one: showing most changing variable in detail
            %
            if isnan(Opt.livePlotFig) || ~Opt.bifurcation.trace
                hs3 = subplot(2,3,4:5);
                plDet = plot(Path.lAll, Path.varAll(interesting,:),'LineWidth', 2);
                set(plDet, {'Color'}, colors(interesting));
                grid on;
                msg = ['interesting variable: $v_{',num2str(interesting),'}$'];
                if isnan(Opt.plotVarOfInterest) msg = [msg, ' (most changing)']; end
                title(msg,'interpreter','latex');
                xlabel('$\lambda$','interpreter','latex');
                ylabel(['$v_{',num2str(interesting),'}$'],'interpreter','latex');
                xlim([Path.lAll-2*dsim1, Path.lAll+2*dsim1]);
                vDiff = abs(vPre(interesting,end) - vPre(interesting,1));
                if 1.1*vDiff > 2*dsim1
                    ylim([Path.varAll(interesting)-vDiff, Path.varAll(interesting)+vDiff]);
                else
                    ylim([Path.varAll(interesting)-2*dsim1, Path.varAll(interesting)+2*dsim1]);
                end
                hold on;
                plCorAssist = plot(lCorAssist(interesting,:), vCorAssist(interesting,:), 'r--', 'LineWidth', 1);
                plCor = plot(lCor(interesting,:), vCor(interesting,:), 'r', 'LineWidth', 1);
                plPre = plot(lPre(interesting,:), vPre(interesting,:), 'k-.', 'LineWidth', 1);
%                 hold off;
                set(hs3,'Tag','lowerleft');
            else
                hs3 = findobj('Tag', 'lowerleft');
                set(hs3, 'NextPlot', 'add');
                plDet = plot(hs3, Path.lAll, Path.varAll(interesting,:),'LineWidth', 2);
                set(plDet, {'Color'}, colors(interesting));
            end
            %
            % second one: showing Path.sAll over Counter.loop
            %
            if isnan(Opt.livePlotFig) || ~Opt.bifurcation.trace
                hs4 = subplot(2,3,6);
                hold off;
                plS = plot(Counter.loop, Path.sAll,'LineWidth', 2, 'Color', 'r');
                grid on;
                title('arc length $s$','interpreter','latex');
                xlabel('loop counter','interpreter','latex');
                ylabel('$s_{\mathrm{all}}$','interpreter','latex');
                set(hs4,'Tag','lowerright');
            else
                hs4 = findobj('Tag', 'lowerright');
                set(hs4, 'NextPlot', 'add');
                plS = plot(hs4, Counter.loop, Path.sAll,'LineWidth', 2, 'Color', 'r');
            end
            %
            Plot.fig = fig; Plot.pl = pl; Plot.plIt = plIt;
            Plot.plDet = plDet; Plot.plS = plS; Plot.plCorAssist = plCorAssist;
            Plot.plCor = plCor; Plot.plPre = plPre;
            Plot.plCurr = plCurr;
            if isnan(Opt.livePlotFig)
                Opt.livePlotFig = fig.Number; % reference to existing fig
            end
        elseif Bifurcation.flag == -1
            %% final change in live plot
            %
            set(0, 'currentfigure', Plot.fig);
            %
            %% calc limits
            dl = abs(Path.lAll(end)-Path.lAll(end-1));
            dv = abs(Path.varAll(interesting,end)-Path.varAll(interesting,end-1));
            %
            %% upper subplots
            %
            % first one
            %
            subplot(2,3,1:2);
            %
            %% save plot data
            %
            newXData = cell(length(Opt.plotVarsIndex), 1);
            newYData = cell(length(Opt.plotVarsIndex), 1);
%             newXDataCurr = cell(length(Opt.plotVarsIndex), 1);
%             newYDataCurr = cell(length(Opt.plotVarsIndex), 1);
            for k = 1:length(Opt.plotVarsIndex)
                newXData{k,1} = Path.lAll;
                newYData{k,1} = Path.varAll(Opt.plotVarsIndex(k),:);
%                 newXDataCurr{k,1} = Path.lAll(end);
%                 newYDataCurr{k,1} = Path.varAll(Opt.plotVarsIndex(k),end);
            end
            %
            %% plot new plot Data
            %
            set(Plot.pl, {'XData'}, newXData, {'YData'},  newYData);
%             set(Plot.plCurr, {'XData'}, newXDataCurr, {'YData'},  newYDataCurr);
            %
            %% add third information to plot (current step)
            for k = 1:length(Opt.plotVarsIndex)
                row = dataTipTextRow('Step',0:(length(Path.lAll)-1),'%d');
                Plot.pl(k).DataTipTemplate.DataTipRows(end+1) = row;
            end
            %
            %% adjust x axis
            if Opt.bifurcation.trace
                axis([Info.lStart, Info.lEnd, min(min(Path.varAll)), max(max(Path.varAll))]);
            else
                xlim([max([lLu(1),lMax(1)-dl0]),min([lLu(2),lMax(2)+dl0])]);
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
            Plot.plDet.XData = Path.lAll;
            Plot.plDet.YData = Path.varAll(interesting,:);
            set(Plot.plDet, {'Color'}, {plot.getRGB(interesting,numPl,1)});
            msg = ['interesting variable: $v_{',num2str(interesting),'}$'];
            if isnan(Opt.plotVarOfInterest) msg = [msg, ' (most changing)']; end
            title(msg,'interpreter','latex');
            ylabel(['$v_{',num2str(interesting),'}$'],'interpreter','latex');
            %
            xlim([Path.lAll(end-1)-2*dl, Path.lAll(end-1)+2*dl]);
            vDiff = abs(vPre(interesting,end) - vPre(interesting,1));
            if 1.1*vDiff > 2*dv
                ylim([Path.varAll(interesting,end-1)-vDiff, Path.varAll(interesting,end-1)+vDiff]);
            else
                ylim([Path.varAll(interesting,end-1)-2*dv, Path.varAll(interesting,end-1)+2*dv]);
            end
            Plot.plCorAssist.XData = lCorAssist(interesting,:);
            Plot.plCorAssist.YData = vCorAssist(interesting,:);
            Plot.plCor.XData = lCor(interesting,:);
            Plot.plCor.YData = vCor(interesting,:);
            Plot.plPre.XData = lPre(interesting,:);
            Plot.plPre.YData = vPre(interesting,:);
            %
            % second one: Path.sAll over Counter.loop
            %
                % no update needed
            %
        else
            set(0, 'currentfigure', Plot.fig);
            %% calc limits
            dl = abs(Path.lAll(end)-Path.lAll(end-1));
            dv = abs(Path.varAll(interesting,end)-Path.varAll(interesting,end-1));
            %
            %% creater upper subplots
            %
            % first one
            %
            subplot(2,3,1:2);
            %% save plot data
            %
            newXData = cell(length(Opt.plotVarsIndex), 1);
            newYData = cell(length(Opt.plotVarsIndex), 1);
            newXDataCurr = cell(length(Opt.plotVarsIndex), 1);
            newYDataCurr = cell(length(Opt.plotVarsIndex), 1);
            for k = 1:length(Opt.plotVarsIndex)
                newXData{k,1} = Path.lAll;
                newYData{k,1} = Path.varAll(Opt.plotVarsIndex(k),:);
                newXDataCurr{k,1} = Path.lAll(end);
                newYDataCurr{k,1} = Path.varAll(Opt.plotVarsIndex(k),end);
            end
            %
            %% plot new plot Data
            %
            set(Plot.pl, {'XData'}, newXData, {'YData'},  newYData);
            set(Plot.plCurr, {'XData'}, newXDataCurr, {'YData'},  newYDataCurr);
            %
            %% mark bifurcation points
            %
            if aux.ison(Opt.bifurcation) && Bifurcation.flag 
               if ~isempty(Bifurcation.bif)
                   hold on;
                   if Bifurcation.bif(2,end) == 0 % brach point Bifurcation.bif
                       plot(Path.lAll(Bifurcation.bif(1,end)),Path.varAll(Opt.plotVarsIndex,Bifurcation.bif(1,end)),'rx','LineWidth',2);
                   elseif Bifurcation.bif(2,end) == 1 % fold Bifurcation.bif
                       plot(Path.lAll(Bifurcation.bif(1,end)),Path.varAll(Opt.plotVarsIndex,Bifurcation.bif(1,end)),'ro','LineWidth',2);
                   elseif isnan(Bifurcation.bif(2,end)) % Bifurcation.bif, but no further information
                       plot(Path.lAll(Bifurcation.bif(1,end)),Path.varAll(Opt.plotVarsIndex,Bifurcation.bif(1,end)),'rs','LineWidth',2);
                   end
                   hold off;
               end
            end
            %
            %% adjust x axis
            %
            xlim([max([lLu(1),lMax(1)-dl0]),min([lLu(2),lMax(2)+dl0])]);
            %
            %% creat lower subplots
            %
            % first one: showing iteration steps
            %
            subplot 233
            %% save plot data
            %
            Plot.plIt.XData = [Plot.plIt.XData, Counter.loop];
            Plot.plIt.YData = [Plot.plIt.YData, iterations];
            xl = xlim();
            xlim([xl(1), Counter.loop]);
            %
            %
            % second one: showing most changing variable in detail
            %
            subplot(2,3,4:5);
            Plot.plDet.XData = Path.lAll;
            Plot.plDet.YData = Path.varAll(interesting,:);
            set(Plot.plDet, {'Color'}, {plot.getRGB(interesting,numPl,1)});
            msg = ['interesting variable: $v_{',num2str(interesting),'}$'];
            if isnan(Opt.plotVarOfInterest) msg = [msg, ' (most changing)']; end
            title(msg,'interpreter','latex');
            ylabel(['$v_{',num2str(interesting),'}$'],'interpreter','latex');
            %
            xlim([Path.lAll(end-1)-2*dl, Path.lAll(end-1)+2*dl]);
            vDiff = abs(vPre(interesting,end) - vPre(interesting,1));
            if 1.1*vDiff > 2*dv
                ylim([Path.varAll(interesting,end-1)-vDiff, Path.varAll(interesting,end-1)+vDiff]);
            else
                ylim([Path.varAll(interesting,end-1)-2*dv, Path.varAll(interesting,end-1)+2*dv]);
            end
            Plot.plCorAssist.XData = lCorAssist(interesting,:);
            Plot.plCorAssist.YData = vCorAssist(interesting,:);
            Plot.plCor.XData = lCor(interesting,:);
            Plot.plCor.YData = vCor(interesting,:);
            Plot.plPre.XData = lPre(interesting,:);
            Plot.plPre.YData = vPre(interesting,:);
            %
            %
            % third one: Path.sAll over Counter.loop
            %
            subplot 236
            Plot.plS.XData = [Plot.plS.XData, Counter.loop];
            Plot.plS.YData = [Plot.plS.YData, Path.sAll(end)];
            xlim([xl(1), Counter.loop]);
            %
            %            
        end
    elseif Opt.plot.threeDim
        %% plot threeDim
        if numPl+1 > 3
            aux.printLine(Opt,'--> 3D plot only works with two variables.\nFirst two are selected! Consider defining plotVarsIndex.\n');
            Opt.plotVarsIndex = [1,2];
            numPl = 1;
        elseif numPl+1 < 3
            error('3D plot only works with at least two variables!');
        end
        if length(Path.lAll) == 1
            if isnan(Opt.livePlotFig) % test for existing figure to plot in
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
                oldFig = false;
            elseif ~ishandle(Opt.livePlotFig)
                error('No such figure exists!');
            else
                fig = figure(Opt.livePlotFig); % use existing fig
                hold on; % for new plot
                oldFig = true;
            end
            %% prepare colors       
            color = plot.getRGB(1,numPl,1);
            %% create plot with colors
            pl = plot3(Path.lAll,Path.varAll(Opt.plotVarsIndex(1),:),Path.varAll(Opt.plotVarsIndex(2),:),'-','LineWidth',2);
            set(pl, 'Color', color); hold on;
            plCurr = plot3(Path.lAll,Path.varAll(Opt.plotVarsIndex(1),:),Path.varAll(Opt.plotVarsIndex(2),:),'*','LineWidth',2);
            set(plCurr, 'Color', color); hold off;
            if ~oldFig
                [caz,cel] = view();
                view(caz+180,cel);
            end
            if isnan(Opt.livePlotFig) || ~Opt.bifurcation.trace % test for existing figure to plot in, there must be no new labels or grid
                grid on;
                xlabel('$\lambda$','interpreter','latex');
                ylabel(['$v_{',num2str(Opt.plotVarsIndex(1)),'}$'],'interpreter','latex');
                zlabel(['$v_{',num2str(Opt.plotVarsIndex(2)),'}$'],'interpreter','latex');
                xlim([max([lLu(1),lMax(1)-dl0]),min([lLu(2),lMax(2)+dl0])]);
            end
            Plot.fig = fig;
            Plot.pl = pl;
            Plot.plCurr = plCurr;
            if isnan(Opt.livePlotFig)
                Opt.livePlotFig = fig.Number; % reference to existing fig
            end
        elseif Bifurcation.flag == -1
            %% final change in live plot
            %
            set(0, 'currentfigure', Plot.fig);
            %
            %% save plot data
            %
            newXData = Path.lAll;
            newYData = Path.varAll(Opt.plotVarsIndex(1),:);
            newZData = Path.varAll(Opt.plotVarsIndex(2),:);
%             newXDataCurr = Path.lAll(end);
%             newYDataCurr = Path.varAll(Opt.plotVarsIndex(1),end);
%             newZDataCurr = Path.varAll(Opt.plotVarsIndex(2),end);
            %
            %% plot new plot Data
            %
            set(Plot.pl, 'XData', newXData, 'YData',  newYData, 'ZData', newZData);
%             set(Plot.plCurr, {'XData'}, newXDataCurr, {'YData'},  newYDataCurr);
            %
            %% add third information to plot (current step)
            row = dataTipTextRow('Step',0:(length(Path.lAll)-1),'%d');
            Plot.pl.DataTipTemplate.DataTipRows(end+1) = row;
            %
            %% adjust x axis
            if Opt.bifurcation.trace
                axis([Info.lStart, Info.lEnd, min(min(Path.varAll)), max(max(Path.varAll))]);
            else
                xlim([max([lLu(1),lMax(1)-dl0]),min([lLu(2),lMax(2)+dl0])]);
            end
        else
            set(0, 'currentfigure', Plot.fig);
            %% save plot data
            %
            newXData = Path.lAll;
            newYData = Path.varAll(Opt.plotVarsIndex(1),:);
            newZData = Path.varAll(Opt.plotVarsIndex(2),:);
            newXDataCurr = Path.lAll(end);
            newYDataCurr = Path.varAll(Opt.plotVarsIndex(1),end);
            newZDataCurr = Path.varAll(Opt.plotVarsIndex(2),end);
            %
            %% plot new plot Data
            %
            set(Plot.pl, 'XData', newXData, 'YData',  newYData, 'ZData', newZData);
            set(Plot.plCurr, 'XData', newXDataCurr, 'YData',  newYDataCurr, 'ZData', newZDataCurr);
            %
            %% mark bifurcation points
            %
            if aux.ison(Opt.bifurcation) && Bifurcation.flag 
               if ~isempty(Bifurcation.bif)
                   hold on;
                   if Bifurcation.bif(2,end) == 0 % brach point Bifurcation.bif
                       plot3(Path.lAll(Bifurcation.bif(1,end)),Path.varAll(Opt.plotVarsIndex(1),Bifurcation.bif(1,end)),Path.varAll(Opt.plotVarsIndex(2),Bifurcation.bif(1,end)),'rx','LineWidth',2);
                   elseif Bifurcation.bif(2,end) == 1 % fold Bifurcation.bif
                       plot3(Path.lAll(Bifurcation.bif(1,end)),Path.varAll(Opt.plotVarsIndex(1),Bifurcation.bif(1,end)),Path.varAll(Opt.plotVarsIndex(2),Bifurcation.bif(1,end)),'ro','LineWidth',2);
                   elseif isnan(Bifurcation.bif(2,end)) % Bifurcation.bif, but no further information
                       plot3(Path.lAll(Bifurcation.bif(1,end)),Path.varAll(Opt.plotVarsIndex(1),Bifurcation.bif(1,end)),Path.varAll(Opt.plotVarsIndex(2),Bifurcation.bif(1,end)),'rs','LineWidth',2);
                   end
                   hold off;
               end
            end
            %
            %% adjust x axis
            %
            xlim([max([lLu(1),lMax(1)-dl0]),min([lLu(2),lMax(2)+dl0])]);        
        end
    elseif Opt.plot.dpa
        %% plot dpa
        if length(Path.lAll) == 1
            if isnan(Opt.livePlotFig) % test for existing figure to plot in
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
                oldFig = false;
            elseif ~ishandle(Opt.livePlotFig)
                error('No such figure exists!');
            else
                fig = figure(Opt.livePlotFig); % use existing fig
                hold on; % for new plot
                oldFig = true;
            end
            %% prepare colors
            colors = cell(numPl,1);        
            for k = 1:numPl
                colors(k) = {plot.getRGB(k,numPl,1)};
            end
            %% create plot with colors
            if Opt.dpaGammaVar
                pl = plot3(Path.lAll,Path.varAll(end,:),Path.varAll(1:(end-1),:),'-','LineWidth',2);
                set(pl, {'Color'}, colors); hold on;
                plCurr = plot3(Path.lAll,Path.varAll(end,:),Path.varAll(1:(end-1),:),'*','LineWidth',2);
                set(plCurr,{'Color'},colors);
                plDpa = plot3(NaN,NaN,NaN,'kd','LineWidth',2);
                hold off;
            else
                pl = plot3(Path.lAll,Opt.g0*ones(size(Path.lAll)),Path.varAll,'-','LineWidth',2);
                set(pl, {'Color'}, colors); hold on;
                plCurr = plot3(Path.lAll,Opt.g0*ones(size(Path.lAll)),Path.varAll,'*','LineWidth',2);
                set(plCurr,{'Color'},colors);
                plDpa = plot3(NaN,NaN,NaN,'kd','LineWidth',2);
                hold off;
            end
            if isnan(Opt.livePlotFig) || ~Opt.bifurcation.trace % test for existing figure to plot in, there must be no new labels or grid
                grid on;
                xlabel('$\lambda$','interpreter','latex');
                ylabel('$\gamma$','interpreter','latex');
                zlabel('$v_{i}$','interpreter','latex');
                if Opt.dpaGammaVar
                    ylim([max([lLu(1),lMax(1)-dl0]),min([lLu(2),lMax(2)+dl0])]);
                else
                    xlim([max([lLu(1),lMax(1)-dl0]),min([lLu(2),lMax(2)+dl0])]);
                    ylim(sort(Opt.g0+[0,(Opt.gTarget-Opt.g0)/100]));
                end
            end
            Plot.fig = fig;
            Plot.pl = pl;
            Plot.plCurr = plCurr;
            Plot.plDpa = plDpa;
            if isnan(Opt.livePlotFig)
                Opt.livePlotFig = fig.Number; % reference to existing fig
            end
        elseif Bifurcation.flag == -1
            %% final change in live plot
            %
            set(0, 'currentfigure', Plot.fig);
            %
            %% save plot data
            %
            if Opt.dpaGammaVar
                newXData = Path.varAll(end,:);
                newYData = Path.lAll;
                newZData = Path.varAll(1:(end-1),:);
            else
                if numel(Path.lAll(:,1))==1
                    newXData = Path.lAll;
                    newYData = Opt.g0*ones(size(Path.lAll));
                    newZData = Path.varAll;
                    newXDataDpa = Path.lAll(1,dpaPoints);
                    newYDataDpa = Opt.g0*ones(size(dpaPoints));
                    newZDataDpa = Path.varAll(:,dpaPoints);
                else
                    newXData = Path.lAll(1,:);
                    newYData = Path.lAll(2,:);
                    newZData = Path.varAll;
                    newXDataDpa = Path.lAll(1,dpaPoints);
                    newYDataDpa = Path.lAll(2,dpaPoints);
                    newZDataDpa = Path.varAll(:,dpaPoints);
                end
            end
            %
            %% plot new plot Data
            %
            set(Plot.pl, 'XData', newXData, 'YData',  newYData, {'ZData'}, num2cell(newZData,2));
            if ~Opt.dpaGammaVar
                set(Plot.plDpa, 'XData', newXDataDpa, 'YData',  newYDataDpa, {'ZData'}, num2cell(newZDataDpa,2));
            end
%             set(Plot.plCurr, {'XData'}, newXDataCurr, {'YData'},  newYDataCurr);
            %
            %% add third information to plot (current step)
            row = dataTipTextRow('Step',0:(length(Path.lAll)-1),'%d');
            for ii=1:numel(Plot.pl)
                Plot.pl(ii).DataTipTemplate.DataTipRows(end+1) = row;
            end
            %
            %% adjust axis
            if Opt.dpaGammaVar
                ylim([max([lLu(1),lMax(1)-dl0]),min([lLu(2),lMax(2)+dl0])]);
            else
                xlim([min(newXData),max(newXData)]);
                ylim([Opt.g0,Opt.gTarget]);
            end
        else
            set(0, 'currentfigure', Plot.fig);
            %% save plot data
            %
            if Opt.dpaGammaVar
                newXData = Path.varAll(end,:);
                newYData = Path.lAll;
                newZData = Path.varAll(1:(end-1),:);
                newXDataCurr = Path.varAll(end,end);
                newYDataCurr = Path.lAll(end);
                newZDataCurr = Path.varAll(1:(end-1),end);
                newXDataDpa = [];
                newYDataDpa = [];
                newZDataDpa = [];
            else
                newXData = Path.lAll;
                newYData = Opt.g0*ones(size(Path.lAll));
                newZData = Path.varAll;
                newXDataCurr = Path.lAll(end);
                newYDataCurr = Opt.g0;
                newZDataCurr = Path.varAll(:,end);
                newXDataDpa = Path.lAll(dpaPoints);
                newYDataDpa = Opt.g0*ones(size(dpaPoints));
                newZDataDpa = Path.varAll(:,dpaPoints);
            end
            %
            %% plot new plot Data
            %
            set(Plot.pl, 'XData', newXData, 'YData', newYData, {'ZData'}, num2cell(newZData,2));
            set(Plot.plCurr, 'XData', newXDataCurr, 'YData', newYDataCurr, {'ZData'}, num2cell(newZDataCurr,2));
            if ~Opt.dpaGammaVar && ~isempty(newXDataDpa)
                set(Plot.plDpa, 'XData', newXDataDpa, 'YData', newYDataDpa, {'ZData'}, num2cell(newZDataDpa,2));
            end
            %
            %% mark bifurcation points
            %
            if aux.ison(Opt.bifurcation) && Bifurcation.flag 
               if ~isempty(Bifurcation.bif)
                   hold on;
                   if Opt.dpaGammaVar
                       if Bifurcation.bif(2,end) == 0 % brach point Bifurcation.bif
                           plot3(Path.lAll(Bifurcation.bif(1,end)),Path.varAll(end,Bifurcation.bif(1,end)),Path.varAll(1:(end-1),Bifurcation.bif(1,end)),'rx','LineWidth',2);
                       elseif Bifurcation.bif(2,end) == 1 % fold Bifurcation.bif
                           plot3(Path.lAll(Bifurcation.bif(1,end)),Path.varAll(end,Bifurcation.bif(1,end)),Path.varAll(1:(end-1),Bifurcation.bif(1,end)),'ro','LineWidth',2);
                       elseif isnan(Bifurcation.bif(2,end)) % Bifurcation.bif, but no further information
                           plot3(Path.lAll(Bifurcation.bif(1,end)),Path.varAll(end,Bifurcation.bif(1,end)),Path.varAll(1:(end-1),Bifurcation.bif(1,end)),'rs','LineWidth',2);
                       end
                   else
                       if Bifurcation.bif(2,end) == 0 % brach point Bifurcation.bif
                           plot3(Path.lAll(Bifurcation.bif(1,end)),Opt.g0,Path.varAll(:,Bifurcation.bif(1,end)),'rx','LineWidth',2);
                       elseif Bifurcation.bif(2,end) == 1 % fold Bifurcation.bif
                           plot3(Path.lAll(Bifurcation.bif(1,end)),Opt.g0,Path.varAll(:,Bifurcation.bif(1,end)),'ro','LineWidth',2);
                       elseif isnan(Bifurcation.bif(2,end)) % Bifurcation.bif, but no further information
                           plot3(Path.lAll(Bifurcation.bif(1,end)),Opt.g0,Path.varAll(:,Bifurcation.bif(1,end)),'rs','LineWidth',2);
                       end
                   end
                   hold off;
               end
            end
            %
            %% adjust x axis
            %
            if Opt.dpaGammaVar
                xl = get(Plot.fig.Children,'XLim');
                xlim([min([xl(1),Path.varAll(end,:)]),max([xl(2),Path.varAll(end,:)])]);
                ylim([max([lLu(1),lMax(1)-dl0]),min([lLu(2),lMax(2)+dl0])]);
            else
                xlim([max([lLu(1),lMax(1)-dl0]),min([lLu(2),lMax(2)+dl0])]);
            end
        end
    else
        error('No such plot method!');
    end
    drawnow limitrate;
    %
    %% Stop to show predictor and corrector
    %
    if islogical(Opt.plotPause) && Opt.plotPause
        aux.printLine(Opt,'Press any key to continue or press Ctrl+c to stop...');
        pause
        aux.printLine(Opt,repmat('\b',1,52));
    elseif numel(Path.lAll) >= Opt.plotPause && ~islogical(Opt.plotPause)
        aux.printLine(Opt,'Press any key to continue or press Ctrl+c to stop...');
        pause
        aux.printLine(Opt,repmat('\b',1,52));
    end
    %
end