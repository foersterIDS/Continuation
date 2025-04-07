%% path continuation - plot.livePlot
%
%   Institute of Dynamics and Vibration Research
%   Leibniz University Hannover
%   11.09.2024 - Tido Kubatschek
%
function livePlot(oih,dpaPoints)
    arguments
        oih
        dpaPoints
    end
    %% Evaluation
    % is the normal plot mode used?
    normalPlotFlag = strcmp(oih.opt.plotOptions.mode,'default');
    % flag for initiation
    isInit = false;
    % calculate amount of added or deleted points
    difference = oih.path.nAll - oih.plot.nPlot;
    % Todo: inserted points !!
    if difference > 0
        for kPlot = 1:difference
            %% Calculate axis entries
            [xVar,yVar,zVar] = evaluatePlotFunction(oih,oih.path.nAll-difference+kPlot,normalPlotFlag);
            if ~isempty(zVar)
                dimension = 3;
            else
                dimension = 2;
            end
            if ~isempty(zVar) && range([size(xVar,2),size(yVar,2),size(zVar,2)])>0
                error('xVar, yVar and zVar must have the same number of columns');
            elseif  range([size(xVar,2),size(yVar,2)])>0
                error('xVar and yVar must have the same number of columns');
            end
            % get numPl
            if dimension == 3
                numPl = max([numel(yVar),numel(zVar)]);
            else
                numPl = numel(yVar);
            end
            %% Initialize
            if oih.path.nAll == 1
                useExisting = false;
                if isnan(oih.opt.plotOptions.figure) || ~ishandle(oih.opt.plotOptions.figure) %isnan(oih.opt.plotOptions.figure) % test for existing figure to plot in
                    fig = figure('units', 'normalized', 'position', [0.2,0.3,0.6,0.5]);
                    clf;
                    oih.opt.plotOptions.figure = fig.Number;
                else % use existing fig
                    fig = figure(oih.opt.plotOptions.figure);
                    hold on; % for new plot
                    oih.opt.plotOptions.figure = fig.Number;
                    useExisting = true;
                end
                if oih.opt.plotOptions.isInitial && ~useExisting
                    isInit = true;
                    % create axis and set axis options
                    ax = axes(oih.opt.plotOptions.axisOptions{:}); hold on;
                    % adapt interpreters
                    ax.XAxis.Label.String = oih.opt.plotOptions.XLabel;
                    ax.YAxis.Label.String = oih.opt.plotOptions.YLabel;
                    ax.ZAxis.Label.String = oih.opt.plotOptions.ZLabel;
                    ax.XAxis.Label.Interpreter = "latex";
                    ax.YAxis.Label.Interpreter = "latex";
                    ax.ZAxis.Label.Interpreter = "latex";
                    grid(oih.opt.plotOptions.grid);
                    box on;
                    if dimension == 3
                        [caz,cel] = view(3);
                        if ~strcmp(oih.opt.plotOptions.mode,'dpa')
                            view(caz+180,cel);
                        end
                    end
                end
                % get Colors
                colors = num2cell(plot.getRGB(1:numPl,numPl,1,oih.opt.plotOptions.colormap),2);
                % create plots
                if dimension < 3
                    pl = plot(xVar,yVar,'LineWidth',oih.opt.plotOptions.LineWidth,'LineStyle',oih.opt.plotOptions.LineStyle); hold on;
                    plCurr = plot(xVar,yVar,'o','LineWidth',2,'LineStyle',oih.opt.plotOptions.LineStyle,'MarkerSize',5);
                else
                    pl = plot3(xVar,yVar,zVar,'LineWidth',oih.opt.plotOptions.LineWidth,'LineStyle',oih.opt.plotOptions.LineStyle); hold on;
                    plCurr = plot3(xVar,yVar,zVar,'o','LineWidth',2,'LineStyle',oih.opt.plotOptions.LineStyle,'MarkerSize',5);
                end
                % set colors
                set(pl, {'Color'}, colors);
                set(plCurr, {'Color'}, colors); hold off;
                set(plCurr, {'MarkerFaceColor'}, colors); hold off;
                oih.plot.fig = fig;
                oih.plot.pl = pl;
                oih.plot.plCurr = plCurr;
            else
                %% Normal routine
                set(0, 'currentfigure', oih.plot.fig);
                %% update plots
                if isscalar(xVar)
                    xVarExtended = xVar*ones(numPl,1);
                else
                    xVarExtended = xVar;
                end
                if dimension == 3 && isscalar(yVar)
                    yVarExtended = yVar*ones(numPl,1);
                else
                    yVarExtended = yVar;
                end
                set(oih.plot.pl,{'XData'},num2cell([cell2mat(get(oih.plot.pl,{'XData'})), xVarExtended],2));
                set(oih.plot.pl,{'YData'},num2cell([cell2mat(get(oih.plot.pl,{'YData'})), yVarExtended],2));
                if ~oih.info.finalSolutionPoint
                    set(oih.plot.plCurr,{'XData'},num2cell(xVarExtended,2));
                    set(oih.plot.plCurr,{'YData'},num2cell(yVar,2));
                end
                if dimension == 3
                    set(oih.plot.pl,{'ZData'},num2cell([cell2mat(get(oih.plot.pl,{'ZData'})), zVar],2));
                    if ~oih.info.finalSolutionPoint
                        set(oih.plot.plCurr,{'ZData'},num2cell(zVar,2));
                    end
                end
            end
            %% adjust axis limits
            %
            if normalPlotFlag
                if isInit
                    xl = [min(xVar),max(xVar)];
                    yl = [min(yVar),max(yVar)];
                else
                    xl = xlim();
                    yl = ylim();
                end
                lLu = [min([oih.info.lStart,oih.info.lEnd]),max([oih.info.lStart,oih.info.lEnd])];
                lMax = [min(oih.path.lAll),max(oih.path.lAll)];
                if oih.path.nAll>=2
                    dl0 = abs(max(oih.path.lAll)-min(oih.path.lAll))*0.2;
                else
                    dl0 = abs(oih.info.lEnd-oih.info.lStart);
                end
                newXLimits = [min([xl(1),max([lLu(1),lMax(1)-dl0])]),max([xl(2), min([lLu(2),lMax(2)+dl0])])];
                if strcmp(oih.opt.plotOptions.XScale,'log')
                    newXLimits = max(newXLimits,0);
                end
                xlim(newXLimits);
                % adjust y axis
                yPlotAll = cell2mat(get(oih.plot.pl,{'YData'}));
                yRange = max(range(yPlotAll));
                dY = 0.1*yRange/2;
                newYLimits = [min([yl(1), min(yPlotAll)-dY]), max([yl(2), max(yPlotAll)+dY])];
                if newYLimits(1) == newYLimits(2)
                    deltaY = max(newYLimits(2), 0.01);
                    newYLimits = [newYLimits(1)-deltaY, newYLimits(2)+deltaY];
                end
                if strcmp(oih.opt.plotOptions.YScale,'log')
                    newYLimits = max(newYLimits,0);
                end
                ylim(newYLimits);
            else
                xPlotAll = cell2mat(get(oih.plot.pl,{'XData'}));
                yPlotAll = cell2mat(get(oih.plot.pl,{'YData'}));
                zPlotAll = cell2mat(get(oih.plot.pl,{'ZData'}));
                
                xRange = max(range(xPlotAll));
                yRange = max(range(yPlotAll));
                
                if isInit
                    currentLimits = [min(xVar),max(xVar),min(yVar),max(yVar),min(zVar),max(zVar)];
                else
                    currentLimits = axis();
                end
        
                dX = 0.1*xRange/2;
                dY = 0.1*yRange/2;
        
                newXLimits = [min([currentLimits(1), min(xPlotAll)-dX]), max([currentLimits(2), max(xPlotAll)+dX])];
                newYLimits = [min([currentLimits(3), min(yPlotAll)-dY]), max([currentLimits(4), max(yPlotAll)+dY])];
                if newXLimits(1) == newXLimits(2)
                    deltaX = max(newXLimits(2), 0.01);
                    newXLimits = [newXLimits(1)-deltaX, newXLimits(2)+deltaX];
                end
                if newYLimits(1) == newYLimits(2)
                    deltaY = max(newYLimits(2), 0.01);
                    newYLimits = [newYLimits(1)-deltaY, newYLimits(2)+deltaY];
                end
                if strcmp(oih.opt.plotOptions.XScale,'log')
                    newXLimits = max(newXLimits,0);
                end
                xlim(newXLimits);
        
                if strcmp(oih.opt.plotOptions.mode,'dpa')
                    newYLimits = sort([oih.opt.g0,oih.opt.gTarget]);
                    if strcmp(oih.opt.plotOptions.YScale,'log')
                        newYLimits = abs(newYLimits);
                    end
                    ylim(newYLimits);
                else
                    if strcmp(oih.opt.plotOptions.YScale,'log')
                        newYLimits = max(newYLimits,0);
                    end
                    ylim(newYLimits);
                end
        
                if ~isempty(zPlotAll)
                    zRange = max(range(zPlotAll));
                    dZ = 0.1*zRange/2;
                    newZLimits = [min([currentLimits(5), min(zPlotAll)-dZ]), max([currentLimits(6), max(zPlotAll)+dZ])];
                    if newZLimits(1) == newZLimits(2)
                        deltaZ = max(newZLimits(2), 0.01);
                        newZLimits = [newZLimits(1)-dZ, newZLimits(2)+deltaZ];
                    end
                    if strcmp(oih.opt.plotOptions.ZScale,'log')
                        newZLimits = max(newZLimits,0);
                    end
                    zlim(newZLimits);
                end
            end
            drawnow limitrate;
            if oih.opt.plotOptions.createAnimation && oih.opt.plotOptions.drawFrame()
                oih.opt.plotOptions.writeFrame();
            end
        end
    elseif difference < 0
        % delete points
        toDelete = abs(difference);
        oldX = cell2mat(get(oih.plot.pl,{'XData'}));
        oldY = cell2mat(get(oih.plot.pl,{'YData'}));
        oldZ = cell2mat(get(oih.plot.pl,{'ZData'}));

        set(oih.plot.pl,{'XData'},num2cell(oldX(:,1:(end-toDelete)),2));
        set(oih.plot.pl,{'YData'},num2cell(oldY(:,1:(end-toDelete)),2));
        if ~oih.info.finalSolutionPoint
            set(oih.plot.plCurr,{'XData'},num2cell(oldX(:,(end-toDelete)),2));
            set(oih.plot.plCurr,{'YData'},num2cell(oldY(:,(end-toDelete)),2));
        end

        if ~isempty(oldZ)
            set(oih.plot.pl,{'ZData'},num2cell(oldZ(:,1:(end-toDelete)),2));
            if ~oih.info.finalSolutionPoint
                set(oih.plot.plCurr,{'ZData'},num2cell(oldZ(:,(end-toDelete)),2));
            end
        end
    end
    oih.plot.nPlot = oih.plot.nPlot + difference;
    %% mark bifurcation points
    % First check if bifurcation option is actived
    % Then check if flag is true, meaning there is a new bifurcation point
    % to plot
    if aux.ison(oih.opt.bifurcation) && oih.bifurcation.flag && ~isempty(oih.bifurcation.bif)
        % find index
        bifIdx = oih.bifurcation.bif(1,end);
        % calculate point
        [xVarBif,yVarBif,zVarBif] = evaluatePlotFunction(oih,bifIdx,normalPlotFlag);
        if ~isempty(zVarBif)
            dimension = 3;
        else
            dimension = 2;
        end
        hold on;
        % find corresponding marker
        if oih.bifurcation.bif(2,end) == 0 % branch point oih.bifurcation.bif
            marker = 'x';
        elseif oih.bifurcation.bif(2,end) == 1 % fold oih.bifurcation.bif
            marker = 'o';
        elseif oih.bifurcation.bif(2,end) == 2 % Bifurkation from additional testfunction
            marker = 'd';
        elseif isnan(oih.bifurcation.bif(2,end)) % oih.bifurcation.bif, but no further information
            marker = 's';
        else
            marker = '';
        end
        if ~isempty(marker)
            if dimension == 2
                plBif = plot(xVarBif,yVarBif,'r','LineWidth',2,'Marker',marker);
            else
                plBif = plot3(xVarBif,yVarBif,zVarBif,'r','LineWidth',2,'Marker',marker);
            end
        else
            plBif = [];
        end
        hold off;
        if ~isempty(plBif)
            nBifPlot = numel(oih.plot.plBif);
            if nBifPlot == 0
                oih.plot.plBif{1} = plBif;
            else
                oih.plot.plBif{nBifPlot+1} = plBif;
            end
        end
    end
    %% mark dpa points
    % First check if dpa option is actived
    % Then check if dpaPoints is not empty, meaning there are points to
    % plot
    if aux.ison(oih.opt.dpa) && ~isempty(dpaPoints)
        % find index
        dpaIdx = dpaPoints(end);
        % calculate point
        [xVarBif,yVarBif,zVarBif] = evaluatePlotFunction(oih,dpaIdx,normalPlotFlag);
        if ~isempty(zVarBif)
            dimension = 3;
        else
            dimension = 2;
        end
        hold on;
        if dimension == 2
            plBif = plot(xVarBif,yVarBif,'r','LineWidth',2,'Marker','o');
        else
            plBif = plot3(xVarBif,yVarBif,zVarBif,'r','LineWidth',2,'Marker','o');
        end
        hold off;
        nBifPlot = numel(oih.plot.plBif);
        if nBifPlot == 0
            oih.plot.plBif{1} = plBif;
        else
            oih.plot.plBif{nBifPlot+1} = plBif;
        end
    end
    %
    %% Stop to show if wanted
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
end

function [xVar,yVar,zVar] = evaluatePlotFunction(oih,evaluationIdx,normalPlotFlag)
    if strcmp(oih.opt.plotOptions.mode,'dpa')
        if oih.opt.dpaGammaVar
            xVar = oih.path.varAll(end,evaluationIdx);
            yVar = oih.path.lAll(evaluationIdx);
            zVar = oih.opt.plotOptions.dpaYFunction(oih.path.varAll(1:(end-1),evaluationIdx));
        else
            if isscalar(oih.path.lAll(:,1))
                xVar = oih.path.lAll(evaluationIdx);
                yVar = oih.opt.g0;
                zVar = oih.opt.plotOptions.dpaYFunction(oih.path.varAll(:,evaluationIdx));
            else
                xVar = oih.path.lAll(1,evaluationIdx);
                yVar = oih.path.lAll(2,evaluationIdx);
                zVar = oih.opt.plotOptions.dpaYFunction(oih.path.varAll(:,evaluationIdx));
            end
        end
    else
        if oih.opt.plotOptions.usePath
            [xVar,yVar,zVar] = oih.opt.plotOptions.calcXYZ([],[],'Path',oih.path);
        else
            [xVar,yVar,zVar] = oih.opt.plotOptions.calcXYZ(oih.path.lAll(evaluationIdx), oih.path.varAll(:,evaluationIdx));
        end
    end
    if normalPlotFlag
        if ~isempty(oih.opt.plotOptions.plotVarsIndex)
            yVar = yVar(oih.opt.plotOptions.plotVarsIndex);
        end
    end
    % correct in case of log scale
    if strcmp(oih.opt.plotOptions.XScale,'log')
        xVar = abs(xVar);
    end
    if strcmp(oih.opt.plotOptions.YScale,'log')
        yVar = abs(yVar);
    end
    if strcmp(oih.opt.plotOptions.ZScale,'log')
        zVar = abs(zVar);
    end
end