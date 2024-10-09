%% path continuation - plot.livePlot
%
%   Institute of Dynamics and Vibration Research
%   Leibniz University Hannover
%   30.09.2020 - Tido Kubatschek
%   07.08.2022 - Anna Lefken
%
function livePlot(oih, ds, dsim1, iterations, funPredictor, sPredictor, dpaPoints)
    lLu = [min([oih.info.lStart,oih.info.lEnd]),max([oih.info.lStart,oih.info.lEnd])];
    lMax = [min(oih.path.lAll),max(oih.path.lAll)];
    if oih.path.nAll>=2
        dl0 = abs(max(oih.path.lAll)-min(oih.path.lAll))*0.2;
    else
        dl0 = abs(oih.info.lEnd-oih.info.lStart);
    end
    numPl = numel(oih.opt.plotVarsIndex);
    if oih.opt.dpaGammaVar
        numPl = numPl-1;
    end
    if nargin<12
        dpaPoints = [];
    end
%     if oih.opt.bifurcation.trace
%         oih.opt.livePlotFig = NaN;
%     end
    if (oih.opt.plot.basic || oih.opt.plot.semilogx || oih.opt.plot.semilogy || oih.opt.plot.loglog)
        if isscalar(oih.path.lAll)
            
            if isnan(oih.opt.livePlotFig) % test for existing figure to plot in
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
                fig = figure(oih.opt.livePlotFig); % use existing fig
                hold on; % for new plot
            end
            %% prepare colors
            colors = cell(numPl,1);        
            for k = 1:numPl
                colors(k) = {plot.getRGB(k,numPl,1)};
            end
            %% create plot with colors
            pl = plot(oih.path.lAll,oih.path.varAll(oih.opt.plotVarsIndex,:),'-','LineWidth',2);
            set(pl, {'Color'}, colors); hold on;
            plCurr = plot(oih.path.lAll,oih.path.varAll(oih.opt.plotVarsIndex,:),'*','LineWidth',2);
            set(plCurr, {'Color'}, colors); hold off;
            if oih.opt.plot.semilogx || oih.opt.plot.loglog
                set(gca, 'XScale', 'log')
                % check whether all lambda are positive
                if sum(oih.path.lAll<0)
                    error('To use semilogx-/loglog-scale lamdba must not contain negative values!');
                end
            end
            if oih.opt.plot.semilogy || oih.opt.plot.loglog
                set(gca, 'YScale', 'log')
                % check whether all variables are positive
                if sum(oih.path.varAll<0)
                    error('To use semilogy-/loglog-scale variables must not contain any negative values! Consinder using the abs()-function.')
                end
            end           
            
            if isnan(oih.opt.livePlotFig) || ~oih.opt.bifurcation.trace % test for existing figure to plot in, there must be no new labels or grid
                grid on;
                xlabel('$\lambda$','interpreter','latex');
                ylabel('$v_{i}$','interpreter','latex');
                xlim([max([lLu(1),lMax(1)-dl0]),min([lLu(2),lMax(2)+dl0])]);
            end
            oih.plot.fig = fig;
            oih.plot.pl = pl;
            oih.plot.plCurr = plCurr;
            if isnan(oih.opt.livePlotFig)
                oih.opt.livePlotFig = fig.Number; % reference to existing fig
                oih.optIsSet.livePlotFig = true;
            end
            
        elseif oih.info.finalSolutionPoint
            %% final change in live plot
            %
            set(0, 'currentfigure', oih.plot.fig);
            %
            %% save plot data
            %
            newXData = cell(length(oih.opt.plotVarsIndex), 1);
            newYData = cell(length(oih.opt.plotVarsIndex), 1);
            newXDataCurr = cell(length(oih.opt.plotVarsIndex), 1);
            newYDataCurr = cell(length(oih.opt.plotVarsIndex), 1);
            for k = 1:length(oih.opt.plotVarsIndex)
                newXData{k,1} = oih.path.lAll;
                newYData{k,1} = oih.path.varAll(oih.opt.plotVarsIndex(k),:);
%                 newXDataCurr{k,1} = oih.path.lAll(end);
%                 newYDataCurr{k,1} = oih.path.varAll(oih.opt.plotVarsIndex(k),end);
            end
            %
            %% plot new plot Data
            %
            set(oih.plot.pl, {'XData'}, newXData, {'YData'},  newYData);
%             set(oih.plot.plCurr, {'XData'}, newXDataCurr, {'YData'},  newYDataCurr);
            %
            %% add third information to plot (current step)
            for k = 1:length(oih.opt.plotVarsIndex)
                row = dataTipTextRow('Step',0:(oih.path.nAll-1),'%d');
                oih.plot.pl(k).DataTipTemplate.DataTipRows(end+1) = row;
            end
            %
            %% adjust x axis
            if oih.opt.bifurcation.trace
                curLimits = axis();
                axis([min([curLimits(1),oih.info.lStart]), max([curLimits(2),oih.info.lEnd]), min([curLimits(3),min(oih.path.varAll)]), max([curLimits(4),max(oih.path.varAll)])]);
            else
                xlim([max([lLu(1),lMax(1)-dl0]),min([lLu(2),lMax(2)+dl0])]);
            end
        else
            set(0, 'currentfigure', oih.plot.fig);
            %% save plot data
            %
            newXData = cell(length(oih.opt.plotVarsIndex), 1);
            newYData = cell(length(oih.opt.plotVarsIndex), 1);
            newXDataCurr = cell(length(oih.opt.plotVarsIndex), 1);
            newYDataCurr = cell(length(oih.opt.plotVarsIndex), 1);
            for k = 1:length(oih.opt.plotVarsIndex)
                newXData{k,1} = oih.path.lAll;
                newYData{k,1} = oih.path.varAll(oih.opt.plotVarsIndex(k),:);
                newXDataCurr{k,1} = oih.path.lAll(end);
                newYDataCurr{k,1} = oih.path.varAll(oih.opt.plotVarsIndex(k),end);
            end
            %
            %% plot new plot Data
            %
            set(oih.plot.pl, {'XData'}, newXData, {'YData'},  newYData);
            set(oih.plot.plCurr, {'XData'}, newXDataCurr, {'YData'},  newYDataCurr);
            %
            %% mark bifurcation points
            %
            if aux.ison(oih.opt.bifurcation) && oih.bifurcation.flag 
               if ~isempty(oih.bifurcation.bif)
                   hold on;
                   if oih.bifurcation.bif(2,end) == 0 % brach point oih.bifurcation.bif
                       plot(oih.path.lAll(oih.bifurcation.bif(1,end)),oih.path.varAll(oih.opt.plotVarsIndex,oih.bifurcation.bif(1,end)),'rx','LineWidth',2);
                   elseif oih.bifurcation.bif(2,end) == 1 % fold oih.bifurcation.bif
                       plot(oih.path.lAll(oih.bifurcation.bif(1,end)),oih.path.varAll(oih.opt.plotVarsIndex,oih.bifurcation.bif(1,end)),'ro','LineWidth',2);
                   elseif oih.bifurcation.bif(2,end) == 2 % Bifurkation from additional testfunction
                       plot(oih.path.lAll(oih.bifurcation.bif(1,end)),oih.path.varAll(oih.opt.plotVarsIndex,oih.bifurcation.bif(1,end)),'rd','LineWidth',2);
                   elseif isnan(oih.bifurcation.bif(2,end)) % oih.bifurcation.bif, but no further information
                       plot(oih.path.lAll(oih.bifurcation.bif(1,end)),oih.path.varAll(oih.opt.plotVarsIndex,oih.bifurcation.bif(1,end)),'rs','LineWidth',2);
                   end
                   hold off;
               end
            end
            %
            %% adjust x axis
            %
            if oih.opt.bifurcation.trace
                curLimits = axis();
                axis([min([curLimits(1),oih.info.lStart]), max([curLimits(2),oih.info.lEnd]), min([curLimits(3),min(oih.path.varAll)]), max([curLimits(4),max(oih.path.varAll)])]);
            else
                xlim([max([lLu(1),lMax(1)-dl0]),min([lLu(2),lMax(2)+dl0])]);
            end
        end
%     elseif oih.opt.plot.detail
%         %% find interesting plot variable
%         if isnan(oih.opt.plotVarOfInterest)
%             %% find most changing
%             nChange = 20;
%             [~,interesting] = max(mean(abs(diff(oih.path.varAll(oih.opt.plotVarsIndex,:),1,2))')');
%         else
%             %% take input by user
%             interesting = oih.opt.plotVarOfInterest;
%         end
%         %
%         %% get corrector
%         [lco, vco] = plot.drawCorrector(oih, dsim1);
%         lCorAssist = lco{1};
%         lCor = lco{2};
%         vCorAssist = vco{1};
%         vCor = vco{2};
%         %
%         %% get predictor
%         if nargin>11
%             sPre = linspace(0,sPredictor,50);
%             xPre = funPredictor(sPre);
%             lPre = kron(ones(oih.info.nv,1),xPre(end,:));
%             vPre = xPre(1:(end-1),:);
%         else
%             lPre = NaN(oih.info.nv,1);
%             vPre = NaN(oih.info.nv,1);
%         end
%         %
%         if isscalar(oih.path.lAll)
%             if isnan(oih.opt.livePlotFig) % test for existing figure to plot in
%                 fig = figure('units', 'normalized', 'position', [0.05,0.1,0.9,0.8]); % create new fig
%                 clf;
%             else
%                 fig = figure(oih.opt.livePlotFig); % use existing fig
%                 hold on; % for new plot
%             end
% 
%             colors = cell(numPl,1);        
%             for k = 1:numPl
%                 colors(k) = {plot.getRGB(k,numPl,1)};
%             end
%             %% creater upper subplot
%             %
%             if isnan(oih.opt.livePlotFig) || ~oih.opt.bifurcation.trace
%                 hs1 = subplot(2,3,1:2);
%                 hold on;
%                 pl = plot(oih.path.lAll,oih.path.varAll(oih.opt.plotVarsIndex,:),'-','LineWidth',2);
%                 set(pl, {'Color'}, colors);
%                 plCurr = plot(oih.path.lAll,oih.path.varAll(oih.opt.plotVarsIndex,:),'*','LineWidth',2);
%                 set(plCurr, {'Color'}, colors);
%                 grid on;
%                 title('path continuation','interpreter','latex');
%                 xlabel('$\lambda$','interpreter','latex');
%                 ylabel('$v_{i}$','interpreter','latex');
%                 xlim([max([lLu(1),lMax(1)-dl0]),min([lLu(2),lMax(2)+dl0])]);
%                 set(hs1,'Tag','upperleft');
%             else
%                 hs1 = findobj('Tag', 'upperleft');
%                 set(hs1, 'NextPlot', 'add');
%                 pl = plot(hs1, oih.path.lAll,oih.path.varAll(oih.opt.plotVarsIndex,:),'-','LineWidth',2);
%                 set(pl, {'Color'}, colors);
%                 plCurr = plot(oih.path.lAll,oih.path.varAll(oih.opt.plotVarsIndex,:),'*','LineWidth',2);
%                 set(plCurr, {'Color'}, colors);
%             end
%             %
%             % second one: showing iteration steps
%             %
%             if isnan(oih.opt.livePlotFig) || ~oih.opt.bifurcation.trace
%                 hs2 = subplot(2,3,3);
%                 plIt = plot(oih.counter.loop, iterations, 'r-x', 'LineWidth', 2);
%                 grid on;
%                 title('needed iterations per loop','interpreter','latex');
%                 xlabel('loop counter','interpreter','latex');
%                 ylabel('iterations','interpreter','latex');
%                 set(hs2,'Tag','upperright');
%             else
%                 hs2 = findobj('Tag', 'upperright');
%                 set(hs2, 'NextPlot', 'replacechildren'); %% maybe renew ???
%                 plIt = plot(hs2, oih.counter.loop, iterations, 'r-x', 'LineWidth', 2);
%             end
%             %
%             %% create lower subplots
%             %
%             % first one: showing most changing variable in detail
%             %
%             if isnan(oih.opt.livePlotFig) || ~oih.opt.bifurcation.trace
%                 hs3 = subplot(2,3,4:5);
%                 plDet = plot(oih.path.lAll, oih.path.varAll(interesting,:),'LineWidth', 2);
%                 set(plDet, {'Color'}, colors(interesting));
%                 grid on;
%                 msg = ['interesting variable: $v_{',num2str(interesting),'}$'];
%                 if isnan(oih.opt.plotVarOfInterest) msg = [msg, ' (most changing)']; end
%                 title(msg,'interpreter','latex');
%                 xlabel('$\lambda$','interpreter','latex');
%                 ylabel(['$v_{',num2str(interesting),'}$'],'interpreter','latex');
%                 xlim([oih.path.lAll-2*dsim1, oih.path.lAll+2*dsim1]);
%                 vDiff = abs(vPre(interesting,end) - vPre(interesting,1));
%                 if 1.1*vDiff > 2*dsim1
%                     ylim([oih.path.varAll(interesting)-vDiff, oih.path.varAll(interesting)+vDiff]);
%                 else
%                     ylim([oih.path.varAll(interesting)-2*dsim1, oih.path.varAll(interesting)+2*dsim1]);
%                 end
%                 hold on;
%                 plCorAssist = plot(lCorAssist(interesting,:), vCorAssist(interesting,:), 'r--', 'LineWidth', 1);
%                 plCor = plot(lCor(interesting,:), vCor(interesting,:), 'r', 'LineWidth', 1);
%                 plPre = plot(lPre(interesting,:), vPre(interesting,:), 'k-.', 'LineWidth', 1);
% %                 hold off;
%                 set(hs3,'Tag','lowerleft');
%             else
%                 hs3 = findobj('Tag', 'lowerleft');
%                 set(hs3, 'NextPlot', 'add');
%                 plDet = plot(hs3, oih.path.lAll, oih.path.varAll(interesting,:),'LineWidth', 2);
%                 set(plDet, {'Color'}, colors(interesting));
%             end
%             %
%             % second one: showing oih.path.sAll over oih.counter.loop
%             %
%             if isnan(oih.opt.livePlotFig) || ~oih.opt.bifurcation.trace
%                 hs4 = subplot(2,3,6);
%                 hold off;
%                 plS = plot(oih.counter.loop, oih.path.sAll,'LineWidth', 2, 'Color', 'r');
%                 grid on;
%                 title('arc length $s$','interpreter','latex');
%                 xlabel('loop counter','interpreter','latex');
%                 ylabel('$s_{\mathrm{all}}$','interpreter','latex');
%                 set(hs4,'Tag','lowerright');
%             else
%                 hs4 = findobj('Tag', 'lowerright');
%                 set(hs4, 'NextPlot', 'add');
%                 plS = plot(hs4, oih.counter.loop, oih.path.sAll,'LineWidth', 2, 'Color', 'r');
%             end
%             %
%             oih.plot.fig = fig; oih.plot.pl = pl; oih.plot.plIt = plIt;
%             oih.plot.plDet = plDet; oih.plot.plS = plS; oih.plot.plCorAssist = plCorAssist;
%             oih.plot.plCor = plCor; oih.plot.plPre = plPre;
%             oih.plot.plCurr = plCurr;
%             if isnan(oih.opt.livePlotFig)
%                 oih.opt.livePlotFig = fig.Number; % reference to existing fig
%             end
%         elseif oih.info.finalSolutionPoint
%             %% final change in live plot
%             %
%             set(0, 'currentfigure', oih.plot.fig);
%             %
%             %% calc limits
%             dl = abs(oih.path.lAll(end)-oih.path.lAll(end-1));
%             dv = abs(oih.path.varAll(interesting,end)-oih.path.varAll(interesting,end-1));
%             %
%             %% upper subplots
%             %
%             % first one
%             %
%             subplot(2,3,1:2);
%             %
%             %% save plot data
%             %
%             newXData = cell(length(oih.opt.plotVarsIndex), 1);
%             newYData = cell(length(oih.opt.plotVarsIndex), 1);
% %             newXDataCurr = cell(length(oih.opt.plotVarsIndex), 1);
% %             newYDataCurr = cell(length(oih.opt.plotVarsIndex), 1);
%             for k = 1:length(oih.opt.plotVarsIndex)
%                 newXData{k,1} = oih.path.lAll;
%                 newYData{k,1} = oih.path.varAll(oih.opt.plotVarsIndex(k),:);
% %                 newXDataCurr{k,1} = oih.path.lAll(end);
% %                 newYDataCurr{k,1} = oih.path.varAll(oih.opt.plotVarsIndex(k),end);
%             end
%             %
%             %% plot new plot Data
%             %
%             set(oih.plot.pl, {'XData'}, newXData, {'YData'},  newYData);
% %             set(oih.plot.plCurr, {'XData'}, newXDataCurr, {'YData'},  newYDataCurr);
%             %
%             %% add third information to plot (current step)
%             for k = 1:length(oih.opt.plotVarsIndex)
%                 row = dataTipTextRow('Step',0:(oih.path.nAll-1),'%d');
%                 oih.plot.pl(k).DataTipTemplate.DataTipRows(end+1) = row;
%             end
%             %
%             %% adjust x axis
%             if oih.opt.bifurcation.trace
%                 axis([oih.info.lStart, oih.info.lEnd, min(min(oih.path.varAll)), max(max(oih.path.varAll))]);
%             else
%                 xlim([max([lLu(1),lMax(1)-dl0]),min([lLu(2),lMax(2)+dl0])]);
%             end
%             %
%             % second one: showing iteration steps
%             %
%                % no update needed
%             %
%             %% lower subplots
%             % first one: showing most changing variable in detail
%             %
%             subplot(2,3,4:5)
%             oih.plot.plDet.XData = oih.path.lAll;
%             oih.plot.plDet.YData = oih.path.varAll(interesting,:);
%             set(oih.plot.plDet, {'Color'}, {plot.getRGB(interesting,numPl,1)});
%             msg = ['interesting variable: $v_{',num2str(interesting),'}$'];
%             if isnan(oih.opt.plotVarOfInterest) msg = [msg, ' (most changing)']; end
%             title(msg,'interpreter','latex');
%             ylabel(['$v_{',num2str(interesting),'}$'],'interpreter','latex');
%             %
%             xlim([oih.path.lAll(end-1)-2*dl, oih.path.lAll(end-1)+2*dl]);
%             vDiff = abs(vPre(interesting,end) - vPre(interesting,1));
%             if 1.1*vDiff > 2*dv
%                 ylim([oih.path.varAll(interesting,end-1)-vDiff, oih.path.varAll(interesting,end-1)+vDiff]);
%             else
%                 ylim([oih.path.varAll(interesting,end-1)-2*dv, oih.path.varAll(interesting,end-1)+2*dv]);
%             end
%             oih.plot.plCorAssist.XData = lCorAssist(interesting,:);
%             oih.plot.plCorAssist.YData = vCorAssist(interesting,:);
%             oih.plot.plCor.XData = lCor(interesting,:);
%             oih.plot.plCor.YData = vCor(interesting,:);
%             oih.plot.plPre.XData = lPre(interesting,:);
%             oih.plot.plPre.YData = vPre(interesting,:);
%             %
%             % second one: oih.path.sAll over oih.counter.loop
%             %
%                 % no update needed
%             %
%         else
%             set(0, 'currentfigure', oih.plot.fig);
%             %% calc limits
%             dl = abs(oih.path.lAll(end)-oih.path.lAll(end-1));
%             dv = abs(oih.path.varAll(interesting,end)-oih.path.varAll(interesting,end-1));
%             %
%             %% creater upper subplots
%             %
%             % first one
%             %
%             subplot(2,3,1:2);
%             %% save plot data
%             %
%             newXData = cell(length(oih.opt.plotVarsIndex), 1);
%             newYData = cell(length(oih.opt.plotVarsIndex), 1);
%             newXDataCurr = cell(length(oih.opt.plotVarsIndex), 1);
%             newYDataCurr = cell(length(oih.opt.plotVarsIndex), 1);
%             for k = 1:length(oih.opt.plotVarsIndex)
%                 newXData{k,1} = oih.path.lAll;
%                 newYData{k,1} = oih.path.varAll(oih.opt.plotVarsIndex(k),:);
%                 newXDataCurr{k,1} = oih.path.lAll(end);
%                 newYDataCurr{k,1} = oih.path.varAll(oih.opt.plotVarsIndex(k),end);
%             end
%             %
%             %% plot new plot Data
%             %
%             set(oih.plot.pl, {'XData'}, newXData, {'YData'},  newYData);
%             set(oih.plot.plCurr, {'XData'}, newXDataCurr, {'YData'},  newYDataCurr);
%             %
%             %% mark bifurcation points
%             %
%             if aux.ison(oih.opt.bifurcation) && oih.bifurcation.flag 
%                if ~isempty(oih.bifurcation.bif)
%                    hold on;
%                    if oih.bifurcation.bif(2,end) == 0 % brach point oih.bifurcation.bif
%                        plot(oih.path.lAll(oih.bifurcation.bif(1,end)),oih.path.varAll(oih.opt.plotVarsIndex,oih.bifurcation.bif(1,end)),'rx','LineWidth',2);
%                    elseif oih.bifurcation.bif(2,end) == 1 % fold oih.bifurcation.bif
%                        plot(oih.path.lAll(oih.bifurcation.bif(1,end)),oih.path.varAll(oih.opt.plotVarsIndex,oih.bifurcation.bif(1,end)),'ro','LineWidth',2);
%                    elseif isnan(oih.bifurcation.bif(2,end)) % oih.bifurcation.bif, but no further information
%                        plot(oih.path.lAll(oih.bifurcation.bif(1,end)),oih.path.varAll(oih.opt.plotVarsIndex,oih.bifurcation.bif(1,end)),'rs','LineWidth',2);
%                    end
%                    hold off;
%                end
%             end
%             %
%             %% adjust x axis
%             %
%             xlim([max([lLu(1),lMax(1)-dl0]),min([lLu(2),lMax(2)+dl0])]);
%             %
%             %% creat lower subplots
%             %
%             % first one: showing iteration steps
%             %
%             subplot 233
%             %% save plot data
%             %
%             oih.plot.plIt.XData = [oih.plot.plIt.XData, oih.counter.loop];
%             oih.plot.plIt.YData = [oih.plot.plIt.YData, iterations];
%             xl = xlim();
%             xlim([xl(1), oih.counter.loop]);
%             %
%             %
%             % second one: showing most changing variable in detail
%             %
%             subplot(2,3,4:5);
%             oih.plot.plDet.XData = oih.path.lAll;
%             oih.plot.plDet.YData = oih.path.varAll(interesting,:);
%             set(oih.plot.plDet, {'Color'}, {plot.getRGB(interesting,numPl,1)});
%             msg = ['interesting variable: $v_{',num2str(interesting),'}$'];
%             if isnan(oih.opt.plotVarOfInterest) msg = [msg, ' (most changing)']; end
%             title(msg,'interpreter','latex');
%             ylabel(['$v_{',num2str(interesting),'}$'],'interpreter','latex');
%             %
%             xlim([oih.path.lAll(end-1)-2*dl, oih.path.lAll(end-1)+2*dl]);
%             vDiff = abs(vPre(interesting,end) - vPre(interesting,1));
%             if 1.1*vDiff > 2*dv
%                 ylim([oih.path.varAll(interesting,end-1)-vDiff, oih.path.varAll(interesting,end-1)+vDiff]);
%             else
%                 ylim([oih.path.varAll(interesting,end-1)-2*dv, oih.path.varAll(interesting,end-1)+2*dv]);
%             end
%             oih.plot.plCorAssist.XData = lCorAssist(interesting,:);
%             oih.plot.plCorAssist.YData = vCorAssist(interesting,:);
%             oih.plot.plCor.XData = lCor(interesting,:);
%             oih.plot.plCor.YData = vCor(interesting,:);
%             oih.plot.plPre.XData = lPre(interesting,:);
%             oih.plot.plPre.YData = vPre(interesting,:);
%             %
%             %
%             % third one: oih.path.sAll over oih.counter.loop
%             %
%             subplot 236
%             oih.plot.plS.XData = [oih.plot.plS.XData, oih.counter.loop];
%             oih.plot.plS.YData = [oih.plot.plS.YData, oih.path.sAll(end)];
%             xlim([xl(1), oih.counter.loop]);
%             %
%             %            
%         end
    elseif oih.opt.plot.threeDim
        %% plot threeDim
        if numPl+1 > 3
            aux.printLine(oih,'--> 3D plot only works with two variables.\nFirst two are selected! Consider defining plotVarsIndex.\n');
            oih.opt.plotVarsIndex = [1,2];
            numPl = 1;
        elseif numPl+1 < 3
            error('3D plot only works with at least two variables!');
        end
        if isscalar(oih.path.lAll)
            if isnan(oih.opt.livePlotFig) % test for existing figure to plot in
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
            elseif ~ishandle(oih.opt.livePlotFig)
                error('No such figure exists!');
            else
                fig = figure(oih.opt.livePlotFig); % use existing fig
                hold on; % for new plot
                oldFig = true;
            end
            %% prepare colors       
            color = plot.getRGB(1,numPl,1);
            %% create plot with colors
            pl = plot3(oih.path.lAll,oih.path.varAll(oih.opt.plotVarsIndex(1),:),oih.path.varAll(oih.opt.plotVarsIndex(2),:),'-','LineWidth',2);
            set(pl, 'Color', color); hold on;
            plCurr = plot3(oih.path.lAll,oih.path.varAll(oih.opt.plotVarsIndex(1),:),oih.path.varAll(oih.opt.plotVarsIndex(2),:),'*','LineWidth',2);
            set(plCurr, 'Color', color); hold off;
            if ~oldFig
                [caz,cel] = view();
                view(caz+180,cel);
            end
            if isnan(oih.opt.livePlotFig) || ~oih.opt.bifurcation.trace % test for existing figure to plot in, there must be no new labels or grid
                grid on;
                xlabel('$\lambda$','interpreter','latex');
                ylabel(['$v_{',num2str(oih.opt.plotVarsIndex(1)),'}$'],'interpreter','latex');
                zlabel(['$v_{',num2str(oih.opt.plotVarsIndex(2)),'}$'],'interpreter','latex');
                xlim([max([lLu(1),lMax(1)-dl0]),min([lLu(2),lMax(2)+dl0])]);
            end
            oih.plot.fig = fig;
            oih.plot.pl = pl;
            oih.plot.plCurr = plCurr;
            if isnan(oih.opt.livePlotFig)
                oih.opt.livePlotFig = fig.Number; % reference to existing fig
            end
        elseif oih.info.finalSolutionPoint
            %% final change in live plot
            %
            set(0, 'currentfigure', oih.plot.fig);
            %
            %% save plot data
            %
            newXData = oih.path.lAll;
            newYData = oih.path.varAll(oih.opt.plotVarsIndex(1),:);
            newZData = oih.path.varAll(oih.opt.plotVarsIndex(2),:);
%             newXDataCurr = oih.path.lAll(end);
%             newYDataCurr = oih.path.varAll(oih.opt.plotVarsIndex(1),end);
%             newZDataCurr = oih.path.varAll(oih.opt.plotVarsIndex(2),end);
            %
            %% plot new plot Data
            %
            set(oih.plot.pl, 'XData', newXData, 'YData',  newYData, 'ZData', newZData);
%             set(oih.plot.plCurr, {'XData'}, newXDataCurr, {'YData'},  newYDataCurr);
            %
            %% add third information to plot (current step)
            row = dataTipTextRow('Step',0:(oih.path.nAll-1),'%d');
            oih.plot.pl.DataTipTemplate.DataTipRows(end+1) = row;
            %
            %% adjust x axis
            if oih.opt.bifurcation.trace
                curLimits = axis();
                axis([min([curLimits(1),oih.info.lStart]), max([curLimits(2),oih.info.lEnd]), min([curLimits(3),min(oih.path.varAll)]), max([curLimits(4),max(oih.path.varAll)])]);
            else
                xlim([max([lLu(1),lMax(1)-dl0]),min([lLu(2),lMax(2)+dl0])]);
            end
        else
            set(0, 'currentfigure', oih.plot.fig);
            %% save plot data
            %
            newXData = oih.path.lAll;
            newYData = oih.path.varAll(oih.opt.plotVarsIndex(1),:);
            newZData = oih.path.varAll(oih.opt.plotVarsIndex(2),:);
            newXDataCurr = oih.path.lAll(end);
            newYDataCurr = oih.path.varAll(oih.opt.plotVarsIndex(1),end);
            newZDataCurr = oih.path.varAll(oih.opt.plotVarsIndex(2),end);
            %
            %% plot new plot Data
            %
            set(oih.plot.pl, 'XData', newXData, 'YData',  newYData, 'ZData', newZData);
            set(oih.plot.plCurr, 'XData', newXDataCurr, 'YData',  newYDataCurr, 'ZData', newZDataCurr);
            %
            %% mark bifurcation points
            %
            if aux.ison(oih.opt.bifurcation) && oih.bifurcation.flag 
               if ~isempty(oih.bifurcation.bif)
                   hold on;
                   if oih.bifurcation.bif(2,end) == 0 % brach point oih.bifurcation.bif
                       plot3(oih.path.lAll(oih.bifurcation.bif(1,end)),oih.path.varAll(oih.opt.plotVarsIndex(1),oih.bifurcation.bif(1,end)),oih.path.varAll(oih.opt.plotVarsIndex(2),oih.bifurcation.bif(1,end)),'rx','LineWidth',2);
                   elseif oih.bifurcation.bif(2,end) == 1 % fold oih.bifurcation.bif
                       plot3(oih.path.lAll(oih.bifurcation.bif(1,end)),oih.path.varAll(oih.opt.plotVarsIndex(1),oih.bifurcation.bif(1,end)),oih.path.varAll(oih.opt.plotVarsIndex(2),oih.bifurcation.bif(1,end)),'ro','LineWidth',2);
                   elseif isnan(oih.bifurcation.bif(2,end)) % oih.bifurcation.bif, but no further information
                       plot3(oih.path.lAll(oih.bifurcation.bif(1,end)),oih.path.varAll(oih.opt.plotVarsIndex(1),oih.bifurcation.bif(1,end)),oih.path.varAll(oih.opt.plotVarsIndex(2),oih.bifurcation.bif(1,end)),'rs','LineWidth',2);
                   end
                   hold off;
               end
            end
            %
            %% adjust x axis
            %
            % xlim([max([lLu(1),lMax(1)-dl0]),min([lLu(2),lMax(2)+dl0])]);        
            %% adjust x axis
            if oih.opt.bifurcation.trace
                curLimits = axis();
                axis([min([curLimits(1),oih.info.lStart]), max([curLimits(2),oih.info.lEnd]), min([curLimits(3),min(oih.path.varAll)]), max([curLimits(4),max(oih.path.varAll)])]);
            else
                xlim([max([lLu(1),lMax(1)-dl0]),min([lLu(2),lMax(2)+dl0])]);
            end
        end
    elseif oih.opt.plot.dpa
        %% plot dpa
        if isscalar(oih.path.lAll)
            if isnan(oih.opt.livePlotFig) % test for existing figure to plot in
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
            elseif ~ishandle(oih.opt.livePlotFig)
                error('No such figure exists!');
            else
                fig = figure(oih.opt.livePlotFig); % use existing fig
                hold on; % for new plot
                oldFig = true;
            end
            %% prepare colors
            colors = cell(numPl,1);        
            for k = 1:numPl
                colors(k) = {plot.getRGB(k,numPl,1)};
            end
            %% create plot with colors
            if oih.opt.dpaGammaVar
                pl = plot3(oih.path.lAll,oih.path.varAll(end,:),oih.path.varAll(1:(end-1),:),'-','LineWidth',2);
                set(pl, {'Color'}, colors); hold on;
                plCurr = plot3(oih.path.lAll,oih.path.varAll(end,:),oih.path.varAll(1:(end-1),:),'*','LineWidth',2);
                set(plCurr,{'Color'},colors);
                plDpa = plot3(NaN,NaN,NaN,'kd','LineWidth',2);
                hold off;
            else
                pl = plot3(oih.path.lAll,oih.opt.g0*ones(size(oih.path.lAll)),oih.path.varAll,'-','LineWidth',2);
                set(pl, {'Color'}, colors); hold on;
                plCurr = plot3(oih.path.lAll,oih.opt.g0*ones(size(oih.path.lAll)),oih.path.varAll,'*','LineWidth',2);
                set(plCurr,{'Color'},colors);
                plDpa = plot3(NaN,NaN,NaN,'kd','LineWidth',2);
                hold off;
            end
            if isnan(oih.opt.livePlotFig) || ~oih.opt.bifurcation.trace % test for existing figure to plot in, there must be no new labels or grid
                grid on;
                xlabel('$\lambda$','interpreter','latex');
                ylabel('$\gamma$','interpreter','latex');
                zlabel('$v_{i}$','interpreter','latex');
                if oih.opt.dpaGammaVar
                    ylim([max([lLu(1),lMax(1)-dl0]),min([lLu(2),lMax(2)+dl0])]);
                else
                    xlim([max([lLu(1),lMax(1)-dl0]),min([lLu(2),lMax(2)+dl0])]);
                    ylim(sort(oih.opt.g0+[0,(oih.opt.gTarget-oih.opt.g0)/100]));
                end
            end
            oih.plot.fig = fig;
            oih.plot.pl = pl;
            oih.plot.plCurr = plCurr;
            oih.plot.plDpa = plDpa;
            if isnan(oih.opt.livePlotFig)
                oih.opt.livePlotFig = fig.Number; % reference to existing fig
                oih.optIsSet.livePlotFig = true;
            end
        elseif oih.info.finalSolutionPoint
            %% final change in live plot
            %
            set(0, 'currentfigure', oih.plot.fig);
            %
            %% save plot data
            %
            if oih.opt.dpaGammaVar
                newXData = oih.path.varAll(end,:);
                newYData = oih.path.lAll;
                newZData = oih.path.varAll(1:(end-1),:);
            else
                if isscalar(oih.path.lAll(:,1))
                    newXData = oih.path.lAll;
                    newYData = oih.opt.g0*ones(size(oih.path.lAll));
                    newZData = oih.path.varAll;
                    newXDataDpa = oih.path.lAll(1,dpaPoints);
                    newYDataDpa = oih.opt.g0*ones(size(dpaPoints));
                    newZDataDpa = oih.path.varAll(:,dpaPoints);
                else
                    newXData = oih.path.lAll(1,:);
                    newYData = oih.path.lAll(2,:);
                    newZData = oih.path.varAll;
                    newXDataDpa = oih.path.lAll(1,dpaPoints);
                    newYDataDpa = oih.path.lAll(2,dpaPoints);
                    newZDataDpa = oih.path.varAll(:,dpaPoints);
                end
            end
            %
            %% plot new plot Data
            %
            set(oih.plot.pl, 'XData', newXData, 'YData',  newYData, {'ZData'}, num2cell(newZData,2));
            if ~oih.opt.dpaGammaVar
                set(oih.plot.plDpa, 'XData', newXDataDpa, 'YData',  newYDataDpa, {'ZData'}, num2cell(newZDataDpa,2));
            end
%             set(oih.plot.plCurr, {'XData'}, newXDataCurr, {'YData'},  newYDataCurr);
            %
            %% add third information to plot (current step)
            row = dataTipTextRow('Step',0:(oih.path.nAll-1),'%d');
            for ii=1:numel(oih.plot.pl)
                oih.plot.pl(ii).DataTipTemplate.DataTipRows(end+1) = row;
            end
            %
            %% adjust axis
            if oih.opt.dpaGammaVar
                ylim([max([lLu(1),lMax(1)-dl0]),min([lLu(2),lMax(2)+dl0])]);
            else
                xlim([min(newXData),max(newXData)]);
                ylim([oih.opt.g0,oih.opt.gTarget]);
            end
        else
            set(0, 'currentfigure', oih.plot.fig);
            %% save plot data
            %
            if oih.opt.dpaGammaVar
                newXData = oih.path.varAll(end,:);
                newYData = oih.path.lAll;
                newZData = oih.path.varAll(1:(end-1),:);
                newXDataCurr = oih.path.varAll(end,end);
                newYDataCurr = oih.path.lAll(end);
                newZDataCurr = oih.path.varAll(1:(end-1),end);
                newXDataDpa = [];
                newYDataDpa = [];
                newZDataDpa = [];
            else
                newXData = oih.path.lAll;
                newYData = oih.opt.g0*ones(size(oih.path.lAll));
                newZData = oih.path.varAll;
                newXDataCurr = oih.path.lAll(end);
                newYDataCurr = oih.opt.g0;
                newZDataCurr = oih.path.varAll(:,end);
                newXDataDpa = oih.path.lAll(dpaPoints);
                newYDataDpa = oih.opt.g0*ones(size(dpaPoints));
                newZDataDpa = oih.path.varAll(:,dpaPoints);
            end
            %
            %% plot new plot Data
            %
            set(oih.plot.pl, 'XData', newXData, 'YData', newYData, {'ZData'}, num2cell(newZData,2));
            set(oih.plot.plCurr, 'XData', newXDataCurr, 'YData', newYDataCurr, {'ZData'}, num2cell(newZDataCurr,2));
            if ~oih.opt.dpaGammaVar && ~isempty(newXDataDpa)
                set(oih.plot.plDpa, 'XData', newXDataDpa, 'YData', newYDataDpa, {'ZData'}, num2cell(newZDataDpa,2));
            end
            %
            %% mark bifurcation points
            %
            if aux.ison(oih.opt.bifurcation) && oih.bifurcation.flag 
               if ~isempty(oih.bifurcation.bif)
                   hold on;
                   if oih.opt.dpaGammaVar
                       if oih.bifurcation.bif(2,end) == 0 % brach point oih.bifurcation.bif
                           plot3(oih.path.lAll(oih.bifurcation.bif(1,end)),oih.path.varAll(end,oih.bifurcation.bif(1,end)),oih.path.varAll(1:(end-1),oih.bifurcation.bif(1,end)),'rx','LineWidth',2);
                       elseif oih.bifurcation.bif(2,end) == 1 % fold oih.bifurcation.bif
                           plot3(oih.path.lAll(oih.bifurcation.bif(1,end)),oih.path.varAll(end,oih.bifurcation.bif(1,end)),oih.path.varAll(1:(end-1),oih.bifurcation.bif(1,end)),'ro','LineWidth',2);
                       elseif isnan(oih.bifurcation.bif(2,end)) % oih.bifurcation.bif, but no further information
                           plot3(oih.path.lAll(oih.bifurcation.bif(1,end)),oih.path.varAll(end,oih.bifurcation.bif(1,end)),oih.path.varAll(1:(end-1),oih.bifurcation.bif(1,end)),'rs','LineWidth',2);
                       end
                   else
                       if oih.bifurcation.bif(2,end) == 0 % brach point oih.bifurcation.bif
                           plot3(oih.path.lAll(oih.bifurcation.bif(1,end)),oih.opt.g0,oih.path.varAll(:,oih.bifurcation.bif(1,end)),'rx','LineWidth',2);
                       elseif oih.bifurcation.bif(2,end) == 1 % fold oih.bifurcation.bif
                           plot3(oih.path.lAll(oih.bifurcation.bif(1,end)),oih.opt.g0,oih.path.varAll(:,oih.bifurcation.bif(1,end)),'ro','LineWidth',2);
                       elseif isnan(oih.bifurcation.bif(2,end)) % oih.bifurcation.bif, but no further information
                           plot3(oih.path.lAll(oih.bifurcation.bif(1,end)),oih.opt.g0,oih.path.varAll(:,oih.bifurcation.bif(1,end)),'rs','LineWidth',2);
                       end
                   end
                   hold off;
               end
            end
            %
            %% adjust x axis
            %
            if oih.opt.dpaGammaVar
                xl = get(oih.plot.fig.Children,'XLim');
                xlim([min([xl(1),oih.path.varAll(end,:)]),max([xl(2),oih.path.varAll(end,:)])]);
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
    if islogical(oih.opt.plotPause) && oih.opt.plotPause
        aux.printLine(oih,'Press any key to continue or press Ctrl+c to stop...');
        pause
        aux.printLine(oih,repmat('\b',1,52));
    elseif oih.path.nAll >= oih.opt.plotPause && ~islogical(oih.opt.plotPause)
        aux.printLine(oih,'Press any key to continue or press Ctrl+c to stop...');
        pause
        aux.printLine(oih,repmat('\b',1,52));
    end
    %
end