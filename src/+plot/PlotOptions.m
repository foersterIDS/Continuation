classdef PlotOptions < handle
    %PlotOptions Declares plot options for livePlot feature of continuation
    properties
        plotVarsIndex {mustBeEmptyOrInteger(plotVarsIndex)} = []
        XTicks (1,:) double {mustBeEmptyOrIncreasing(XTicks)}
        YTicks (1,:) double {mustBeEmptyOrIncreasing(YTicks)}
        ZTicks (1,:) double {mustBeEmptyOrIncreasing(ZTicks)}
        XLabel (1,:) char
        YLabel (1,:) char
        ZLabel (1,:) char
        XScale (1,:) char {mustBeMember(XScale,{'linear','log'})}  = 'linear'
        YScale (1,:) char {mustBeMember(YScale,{'linear','log'})}  = 'linear'
        ZScale (1,:) char {mustBeMember(ZScale,{'linear','log'})}  = 'linear'
        figure (1,1) double {mustBeGreaterThanZeroOrNaN(figure,1)} = NaN
        LineStyle (1,:) char {mustBeValidLineStyle(LineStyle)} = '-'
        LineWidth (1,1) double {mustBeGreaterThanZeroOrNaN(LineWidth)} = 2
        colormap (1,:) char = 'default'
        grid (1,:) char = 'on'
        dpaYFunction (1,1) function_handle = @(v) v
    end

    properties (Hidden)
        usePath (1,1) logical = false
    end

    properties (GetAccess=public,SetAccess=private)
        axisOptions
        plotFunction
        mode (1,:) char {mustBeMember(mode,{'default','3d','dpa','custom'})} = 'default'
    end

    properties (Access=private)
        initial = true
    end

    methods
        function obj = PlotOptions(options)
            arguments
                options.plotFunction (1,:) {mustBeEmptyOrFunctionHandle(options.plotFunction)} = []
                options.plotVarsIndex {mustBeEmptyOrInteger(options.plotVarsIndex)} = []
                options.XTicks (1,:) double {mustBeEmptyOrIncreasing(options.XTicks)}
                options.YTicks (1,:) double {mustBeEmptyOrIncreasing(options.YTicks)}
                options.ZTicks (1,:) double {mustBeEmptyOrIncreasing(options.ZTicks)}
                options.XLabel (1,:) char
                options.YLabel (1,:) char
                options.ZLabel (1,:) char
                options.XScale (1,:) char {mustBeMember(options.XScale,{'linear','log'})} = 'linear'
                options.YScale (1,:) char {mustBeMember(options.YScale,{'linear','log'})} = 'linear'
                options.ZScale (1,:) char {mustBeMember(options.ZScale,{'linear','log'})} = 'linear'
                options.idx3d (1,2) double {mustBeInteger,mustBePositive} = [1,2]
                options.mode (1,:) char {mustBeMember(options.mode,{'default','3d','dpa','custom'})} = 'default'
                options.figure (1,1) double {mustBeGreaterThanZeroOrNaN(options.figure,1)} = NaN
                options.LineStyle (1,:) char {mustBeValidLineStyle(options.LineStyle)} = '-'
                options.LineWidth (1,1) double {mustBeGreaterThanZeroOrNaN(options.LineWidth)} = 2
                options.colormap (1,:) char = 'viridis'
                options.grid (1,:) char {mustBeValidGrid(options.grid)} = 'on'
                options.index3D (1,2) double {mustBeInteger} = [1,2]
                options.dpaYFunction (1,1) function_handle = @(v) v
                options.usePath (1,1) logical = false
            end
            %% argument validation
            changeInputOfPlotFunctionToDefault = false;
            usePlotVarsIndex = ~isempty(options.plotVarsIndex);

            if ~isempty(options.plotFunction)
                if usePlotVarsIndex
                    warning('Using plotVarsIndex results in the input of plotFunction to not be used.')
                    changeInputOfPlotFunctionToDefault = true;
                end
                if ~strcmp(options.mode,'custom')
                    if ~strcmp(options.mode,'default')
                        warning('You cannot use mode = %s and define plotFunction. Mode input is not used.',options.mode);
                    end
                    options.mode = 'custom';
                end
            else
                if strcmp(options.mode,{'3d','dpa'}) 
                    if usePlotVarsIndex
                        warning('You cannot use mode = %s and define plotVarsIndex.',options.mode);
                    end
                elseif strcmp(options.mode,'custom')
                    warning('mode input custom ingnored.');
                    options.mode = 'default';
                else
                    changeInputOfPlotFunctionToDefault = true;
                end
            end

            if changeInputOfPlotFunctionToDefault
                options.plotFunction = @(l,v) deal(l,v,[]);
            end

            if strcmp(options.mode,'3d')
                options.plotFunction = @(l,v) plot.plot3d(l,v,'idxY',options.index3D(1),'idxZ',options.index3D(2));
                if ~isfield(options,'XLabel')
                    options.XLabel = '$\lambda$';
                end
                if ~isfield(options,'YLabel')
                    options.YLabel = ['$v_',int2str(options.idx3d(1)),'$'];
                end
                if ~isfield(options,'ZLabel')
                    options.ZLabel = ['$v_',int2str(options.idx3d(2)),'$'];
                end
            elseif strcmp(options.mode,'dpa')
                options.plotFunction = [];
                if ~isfield(options,'XLabel')
                    options.XLabel = '$\lambda$';
                end
                if ~isfield(options,'YLabel')
                    options.YLabel = '$\gamma$';
                end
                if ~isfield(options,'ZLabel')
                    options.ZLabel = '$v_i$';
                end
            else % default
                if ~isfield(options,'XLabel')
                    options.XLabel = '$\lambda$';
                end
                if ~isfield(options,'YLabel')
                    options.YLabel = '$v_i$';
                end
                if ~isfield(options,'ZLabel')
                    options.ZLabel = '';
                end
            end

            %% fill object
            fields = fieldnames(options);
            for kk = 1:numel(fields)
                if isprop(obj,fields{kk})
                    obj.(fields{kk}) = options.(fields{kk});
                end
            end
        end

        %% getter methods
        function axOptions = get.axisOptions(obj)
            axOptions = cell(0,1);
            nOptions = 0;
            possibleFields = {
            'XTicks',...
            'YTicks',...
            'ZTicks',...
            'XScale',...
            'YScale',...
            'ZScale'
            };
            for kk = 1:length(possibleFields)
                if ~isempty(obj.(possibleFields{kk})) && ~any(isnan(obj.(possibleFields{kk})))
                    axOptions{1,nOptions+1} = possibleFields{kk};
                    axOptions{1,nOptions+2} = obj.(possibleFields{kk});
                    nOptions = nOptions + 2;
                end
            end
        end

        function val = isInitial(obj)
            if obj.initial
                val = true;
                obj.initial = false;
            else
                val = false;
            end
        end

        function [x,y,z] = calcXYZ(obj,l,v,NameValueArgs)
            arguments
                obj 
                l 
                v
                NameValueArgs.Path = []
            end
            if obj.usePath
                [x,y,z] = obj.plotFunction(NameValueArgs.Path);
            else
                [x,y,z] = obj.plotFunction(l,v);
            end
        end

        %% Setter Methods
        function obj = changeModeToDPA(obj)
            arguments
                obj (1,1) plot.PlotOptions
            end
            obj.mode = 'dpa';
            obj.plotFunction = [];
        end
    end
end

function mustBeEmptyOrFunctionHandle(var)
    if ~(isa(var,'function_handle') || isempty(var))
        eidType = 'mustBeLogicalOrFunctionHandle:notLogicalOrFunctionHandle';
        msgType = 'Input must be a logical or a function handle defining the axis input.';
        error(eidType,msgType)
    end
end

function mustBeEmptyOrIncreasing(var)
    if ~isempty(var) && (~(all(size(var) == size(unique(var))) && (var == unique(var))))
        eidType = 'mustBeIncreasing:notIncreasing';
        msgType = 'Input must be a real-valued array of increasing values.';
        error(eidType,msgType)
    end
end

function mustBeValidLineStyle(var)
    if ~strcmp(var,{'-','--',':','-.'})
        eidType = 'mustBeValidLineStyle:notAValidLineStyle';
        msgType = 'Input must be a valid linestyle: -, --, :, -. .';
        error(eidType,msgType)
    end
end

function mustBeGreaterThanZeroOrNaN(var,mustBeIntegerFlag)
    arguments
        var 
        mustBeIntegerFlag (1,1) logical = false;
    end
    if mustBeIntegerFlag
        msgType = 'Input must be NaN or a integer value greater than 0.';
    else
        msgType = 'Input must be NaN or a double value greater than 0.';
    end
    if ~isnan(var)
        if var <= 0
            eidType = 'mustBeGreaterThanZeroOrNaN:notZeroOrNaN';
            error(eidType,msgType)
        elseif mustBeIntegerFlag && ~(var == round(var))
            eidType = 'mustBeGreaterThanZeroOrNaN:notZeroOrNaN';
            error(eidType,msgType)
        end
    end
end

function mustBeEmptyOrInteger(var)
    if ~(isempty(var) || all(round(var) == var & var > 0))
        eidType = 'mustBeEmptyOrInteger:notEmptyOrInteger';
        msgType = 'Input must be empty or an array of integers greater than 0.';
        error(eidType,msgType)
    end
end

function mustBeValidGrid(var)
    if ~strcmp(var,{'on','minor','off'})
        eidType = 'mustBeValidGrid:notValidGrid';
        msgType = 'Input must be either: on, off, minor.';
        error(eidType,msgType)
    end
end