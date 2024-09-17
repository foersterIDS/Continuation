%% path continuation - continuation.Path (class)
%
%   Institute of Dynamics and Vibration Research
%   Leibniz University Hannover
%   28.07.2023 - Alwin Förster
%
classdef Path < handle
    %% properties
    properties (GetAccess = public, SetAccess = private)
        bifTestValue            % bifurcation test value
        iterations              % number of solver iterations
        JAll                    % cell array of jacobians
        lAll                    % parameter values (scalar)
        lDir                    % parameter directions (vector)
        lRoot                   % initial parameter
        nL                      % number of parameters
        nVar                    % number of variables
        pathInfoValue           % path info value
        stepBackStatus = false  % status of step back mode
        rateOfContraction       % rate of contraction (solver)
        signDetJRedAll          % sign of determinant of the jacobian (nVar x nVar)
        speedOfContinuation     % speed of continuation
        varAll                  % variable values
        xPredictorAll           % predictor values
    end

    properties (Access = public)
        outputFormat            % output format 'full'/'singleValueForL'
    end

    properties (Access = private)
        JFullTemp               % full jacobian of last evaluation
        oih                     % OptInfoHandle object
        plusStruct = struct('bifTestValue',[],'iterations',[],'J',[],...
                            'pathInfoValue',[],'rateOfContraction',[],...
                            'signDetJ',[],'speedOfContinuation',[],...
                            'var',[],'l',[],'xPredictor',[]);
    end

    properties (Dependent)
        nAll                    % number of points on path
        nX                      % number of vars and ls
        sAll                    % arc length
        xAll                    % variables and parameters
        xPlus                   % variables and parameters of next point in stepback
    end

    %% public methods
    methods
        %% constructor
        function obj = Path(nVar,nL,oih)
            %% arguments
            arguments
                nVar (1,1) double {mustBeInteger,mustBeGreaterThanOrEqual(nVar,1)}
                nL (1,1) double {mustBeInteger,mustBeGreaterThanOrEqual(nL,1)}
                oih (1,1) aux.OptInfoHandle
            end
            %% initialize
            obj.nVar = nVar;
            obj.nL = nL;
            obj.oih = oih;
            if obj.oih.opt.jacobianOut.basic
                obj.JAll = {[],[],[]}; % initial, last, previous
            elseif obj.oih.opt.jacobianOut.full
                obj.JAll = {};
            else
                error('jacobianOut must be full or basic!');
            end
            obj.lAll = zeros(nL,0);
            if oih.info.nl>1
                obj.lRoot = obj.oih.opt.lMult0;
            end
            obj.signDetJRedAll = zeros(1,0);
            obj.varAll = zeros(nVar,0);
            obj.resetOutput();
            %% additional functions
            if obj.oih.optIsSet.pathInfoFunction
                obj.pathInfoValue = zeros(1,0);
            end
            if obj.oih.optIsSet.bifAdditionalTestfunction
                obj.bifTestValue = zeros(1,0);
            end
            %% step size control values
            if obj.oih.stepsizeOptions.iterations
                obj.iterations = zeros(1,0);
            end
            if obj.oih.stepsizeOptions.predictor
                obj.xPredictorAll = zeros(nVar+nL,0);
            end
            if obj.oih.stepsizeOptions.rateOfContraction
                obj.rateOfContraction = zeros(1,0);
            end
            if obj.oih.stepsizeOptions.speedOfContinuation
                obj.speedOfContinuation = zeros(1,0);
            end
        end

        %% getter
        function lAll = get.lAll(obj)
            switch obj.outputFormat
                case 'singleValueForL'
                    if obj.nL==1
                        lAll = obj.lAll;
                    else
                        if isempty(obj.lAll)
                            lAll = zeros(obj.nL,0);
                        else
                            lAll = cumsum([0,sqrt(sum(diff(obj.lAll,1,2).^2,1))]);
                        end
                    end
                case 'full'
                    lAll = obj.lAll;
                otherwise
                    error('Unknown output format.');
            end
        end

        function nAll = get.nAll(obj)
            nAll = numel(obj.lAll(1,:));
        end

        function nX = get.nX(obj)
            nX = obj.nVar+obj.nL;
        end

        function sAll = get.sAll(obj)
            sAll = [0,cumsum(sqrt(sum(diff(obj.xAll,1,2).^2,1)))];
        end

        function xAll = get.xAll(obj)
            xAll = [obj.varAll;obj.lAll];
        end

        function xPlus = get.xPlus(obj)
            xPlus = [obj.plusStruct.var;obj.plusStruct.l];
        end

        %% setter
        function set.outputFormat(obj,oF)
            arguments
                obj (1,1) continuation.Path
                oF (1,:) char {mustBeMember(oF,{'singleValueForL','full'})}
            end
            obj.outputFormat = oF;
        end

        %% general methods
        function addPoint(obj,var,l,J,oih,pos,NameValueArgs)
            %% arguments
            arguments
                obj (1,1) continuation.Path
                var (:,1) double
                l (:,1) double
                J (:,:) double
                oih (1,1) aux.OptInfoHandle
                pos (1,1) double {mustBeInteger,mustBeGreaterThan(pos,0)} = obj.nAll+1
                NameValueArgs.bifTestValue (1,1) double
                NameValueArgs.pathInfoValue (1,1) double
                NameValueArgs.predictor (:,1) double
                NameValueArgs.speedOfContinuation (1,1) double
            end
            validateattributes(var,{'numeric'},{'size',[obj.nVar,1]});
%             validateattributes(l,{'numeric'},{'size',[obj.nL,1]});
            obj.checkInputOptions(NameValueArgs);
            %% transform parameters
            if obj.nL>1
                if pos~=obj.nAll+1
                    error('With nL>1, entries may only be added at the end of the path.');
                end
                if obj.nAll==0
                    %% initial
                    if ~obj.oih.optIsSet.lDirFunction
                        error("'lDirFunction' must be set if more than one parameter is used.");
                    end
                    l = obj.lRoot;
                else
                    %% additional
                    obj.outputFormat = 'full';
                    lLast = obj.lAll(:,pos-1);
                    obj.resetOutput();
                    lDirLast = obj.lDir(:,pos-1);
                    if numel(l)==1
                        %% calc. full l
                        lSVAdd = l;
                        lSVEnd = obj.lAll(:,end);
                        dlSV = lSVAdd-lSVEnd;
                        l = lLast+dlSV*lDirLast;
                        %% calc. full J
                        if size(J,2)<obj.nVar+obj.nL
                            Jred = J;
                            J = obj.JFullTemp;
                            if isempty(J) || ~all(all(Jred(1:obj.nVar,1:obj.nVar)==J(1:obj.nVar,1:obj.nVar)))
                                error('No current jacobian found.')
                            end
                        end
                    else
                        error('Number of parameter input(!) must be 1.');
                    end
                end
            end
            %% add/insert
            obj.varAll = [obj.varAll(:,1:(pos-1)),var,obj.varAll(:,pos:end)];
            obj.outputFormat = 'full';
            obj.lAll = [obj.lAll(:,1:(pos-1)),l,obj.lAll(:,pos:end)];
            obj.resetOutput();
            if obj.nL>1
                if obj.nAll>1
                    ds = norm(diff(obj.xAll(:,end+[-1,0]),1,2));
                else
                    ds = obj.oih.info.ds0;
                end
                lDirTemp = obj.oih.opt.lDirFunction(var,l,J,ds);
                lDirTemp = lDirTemp/norm(lDirTemp);
                obj.lDir = [obj.lDir(:,1:(pos-1)),lDirTemp,obj.lDir(:,pos:end)];
            end
            if obj.oih.opt.jacobianOut.basic
                obj.JAll{3} = obj.JAll{2};
                obj.JAll{2} = J;
            elseif obj.oih.opt.jacobianOut.full
                obj.JAll = [obj.JAll(:,1:(pos-1)),J,obj.JAll(:,pos:end)];
            else
                error('jacobianOut must be full or basic!');
            end
            obj.signDetJRedAll = [obj.signDetJRedAll(:,1:(pos-1)),sign(det(J(1:obj.nVar,1:obj.nVar))),obj.signDetJRedAll(:,pos:end)];
            if obj.oih.optIsSet.bifAdditionalTestfunction
                obj.bifTestValue = [obj.bifTestValue(:,1:(pos-1)),NameValueArgs.bifTestValue,obj.bifTestValue(:,pos:end)];
            end
            if obj.oih.stepsizeOptions.iterations
                obj.iterations = [obj.iterations(:,1:(pos-1)),oih.solver.output.iterations,obj.iterations(:,pos:end)];
            end
            if obj.oih.optIsSet.pathInfoFunction
                obj.pathInfoValue = [obj.pathInfoValue(:,1:(pos-1)),NameValueArgs.pathInfoValue,obj.pathInfoValue(:,pos:end)];
            end
            if obj.oih.stepsizeOptions.predictor
                obj.xPredictorAll = [obj.xPredictorAll(:,1:(pos-1)),NameValueArgs.predictor,obj.xPredictorAll(:,pos:end)];
            end
            if obj.oih.stepsizeOptions.rateOfContraction
                obj.rateOfContraction = [obj.rateOfContraction(:,1:(pos-1)),oih.solver.output.rateOfContraction,obj.rateOfContraction(:,pos:end)];
            end
            if obj.oih.stepsizeOptions.speedOfContinuation
                obj.speedOfContinuation = [obj.speedOfContinuation(:,1:(pos-1)),NameValueArgs.speedOfContinuation,obj.speedOfContinuation(:,pos:end)];
            end
        end

        function addPointAtEnd(obj,var,l,J,oih,NameValueArgs)
            %% arguments
            arguments
                obj (1,1) continuation.Path
                var (:,1) double
                l (:,1) double
                J (:,:) double
                oih (1,1) aux.OptInfoHandle
                NameValueArgs.bifTestValue (1,1) double
                NameValueArgs.pathInfoValue (1,1) double
                NameValueArgs.predictor (:,1) double
                NameValueArgs.speedOfContinuation (1,1) double
            end
            obj.checkInputOptions(NameValueArgs);
            %% pass to addPoint(...)
            pos = obj.nAll+1;
            nva = aux.struct2cellPreserveFieldnames(NameValueArgs);
            obj.addPoint(var,l,J,oih,pos,nva{:});
        end

        function detJv = detJv(obj,NameValueArgs)
            arguments
                obj (1,1) continuation.Path
                NameValueArgs.index (1,1) double {mustBeInteger}
                NameValueArgs.name (1,:) char {mustBeMember(NameValueArgs.name,{'last','previous','initial','plus'})} = 'last'
            end
            if isfield(NameValueArgs,'index')
                Jv = obj.JAll{index}(1:obj.nVar,1:obj.nVar);
            else
                switch NameValueArgs.name
                    case 'last'
                        Jv = obj.JAll{end}(1:obj.nVar,1:obj.nVar);
                    case 'previous'
                        if obj.nAll>1
                            Jv = obj.JAll{end-1}(1:obj.nVar,1:obj.nVar);
                        else
                            Jv = obj.JAll{end}(1:obj.nVar,1:obj.nVar);
                        end
                    case 'initial'
                        Jv = obj.JAll{1}(1:obj.nVar,1:obj.nVar);
                    case 'plus'
                        Jv = obj.plusStruct.J(1:obj.nVar,1:obj.nVar);
                end
            end
            detJv = det(Jv);
        end

        function [varargout] = fOfvarAndlSingleIO(obj,fun,v,l)
            if obj.nL==1
                %% nL = 1
                if obj.oih.opt.jacobian
                    [varargout{1},varargout{2}] = fun(v,l);
                else
                    varargout{1} = fun(v,l);
                end
            else
                %% nL > 1
                % l to lFull
                obj.outputFormat = 'full';
                lLast = obj.lAll(:,end);
                obj.resetOutput();
                lSVAdd = l;
                lSVEnd = obj.lAll(:,end);
                dlSV = lSVAdd-lSVEnd;
                l = lLast+dlSV*obj.lDir(:,end);
                % eval fun
                if obj.oih.opt.jacobian
                    [varargout{1},JFull] = fun(v,l);
                    obj.JFullTemp = JFull;
                    % JFull to J
                    varargout{2} = [JFull(:,1:obj.nVar),JFull(:,obj.nVar+(1:obj.nL))*obj.lDir(:,end)];
                else
                    varargout{1} = fun(v,l);
                end
            end
        end

        function lFull = getFullLfromLSV(obj,lSV)
            if obj.nL==1
                lFull = lSV;
            else
                obj.outputFormat = 'full';
                lLast = obj.lAll(:,end);
                obj.resetOutput();
                lDirLast = obj.lDir(:,end);
                %% calc. full l
                lSVEnd = obj.lAll(:,end);
                if lSV<lSVEnd
                    error('lSV has to be greater than lSVEnd');
                end
                dlSV = lSV-lSVEnd;
                lFull = lLast+dlSV*lDirLast;
            end
        end

        function pathOut = copy(obj)
            pathOut = continuation.Path(obj.nVar,obj.nL,obj.oih);
            props = properties(obj);
            for ii=1:numel(props)
                prop = props{ii};
                try
                    pathOut.(prop) = obj.(prop);
                catch
                    % dependend or private set
                end
            end
        end

        function J = getJacobianByName(obj,jacName,idx1,idx2)
            %% arguments
            arguments
                obj (1,1) continuation.Path
                jacName (1,:) char {mustBeMember(jacName,{'last','previous','initial','plus'})}
                idx1 (1,:) double {mustBeInteger,mustBeGreaterThan(idx1,0)} = 1:obj.nVar
                idx2 (1,:) double {mustBeInteger,mustBeGreaterThan(idx2,0)} = 1:(obj.nVar+obj.nL)
            end
            %% get jacobian
            if obj.oih.opt.jacobianOut.basic
                if numel(obj.JAll)==3
                    switch jacName
                        case 'initial'
                            J = obj.JAll{1};
                        case 'last'
                            J = obj.JAll{2};
                        case 'plus'
                            J = obj.plusStruct.J;
                        case 'previous'
                            if obj.nAll==1
                                J = obj.JAll{2};
                            else
                                J = obj.JAll{3};
                            end
                    end
                else
                    J = [];
                end
            elseif obj.oih.opt.jacobianOut.full
                if obj.nAll>0
                    switch jacName
                        case 'initial'
                            J = obj.JAll{1};
                        case 'last'
                            J = obj.JAll{end};
                        case 'plus'
                            J = obj.plusStruct.J;
                        case 'previous'
                            if obj.nAll==1
                                J = obj.JAll{end};
                            else
                                J = obj.JAll{end-1};
                            end
                    end
                else
                    J = [];
                end
            else
                error('jacobianOut must be full or basic!');
            end
            if ~isinf(idx1)
                idx1(logical((idx1>size(J,1))+(idx1<1))) = [];
                J = J(idx1,:);
            end
            if ~isinf(idx2)
                idx2(logical((idx2>size(J,2))+(idx2<1))) = [];
                J = J(:,idx2);
            end
        end

        function li = lIndex(obj,index)
            arguments
                obj (1,1) continuation.Path
                index (1,:) double {mustBeInteger,mustBeGreaterThanOrEqual(index,1)}
            end
            validateattributes(index,{'double'},'<=',obj.nAll);
            li = obj.lAll(:,index);
        end

        function overwrite(obj,varAll,lAll,bifAll,doClearPath)
            arguments
                obj (1,1) continuation.Path
                varAll (:,:) double
                lAll (:,:) double
                bifAll (:,:) double
                doClearPath (1,1) logical
            end
            if doClearPath
                obj.outputFormat = 'full';
                varAllTemp = obj.varAll;
                lAllTemp = obj.lAll;
                bifTemp = obj.bifTestValue;
                obj.resetOutput();
                obj.clearPath();
            end
            if ~isempty(varAll)
                obj.varAll = varAll;
            elseif doClearPath
                obj.varAll = varAllTemp;
            end
            if ~isempty(lAll)
                obj.lAll = lAll;
            elseif doClearPath
                obj.lAll = lAllTemp;
            end
            if ~isempty(bifAll)
                obj.bifTestValue = bifAll;
            elseif doClearPath
                obj.bifTestValue = bifTemp;
            end
        end

        function remove(obj,idxRmv)
            %% arguments
            arguments
                obj (1,1) continuation.Path
                idxRmv (1,:) double {mustBeInteger,mustBeGreaterThan(idxRmv,0)}
            end
            %% sort idxRmv
            idxRmv = sort(idxRmv);
            %% check nL
            if obj.nL>1 && ~(isempty(idxRmv) || (idxRmv(end)==obj.nAll && all(diff(idxRmv)==1)))
                error('If nL>1, only the last entries in the path may be removed.');
            end
            %% remove
            if obj.oih.opt.jacobianOut.basic
                if sum(idxRmv==obj.nAll)
                    if sum(idxRmv==obj.nAll-1)
                        obj.JAll{2} = [];
                        obj.JAll{3} = [];
                    else
                        obj.JAll{2} = obj.JAll{3};
                        obj.JAll{3} = [];
                    end
                elseif sum(idxRmv==obj.nAll-1)
                    obj.JAll{3} = [];
                end
                if sum(idxRmv==1)
                    obj.JAll{1} = [];
                end
            elseif obj.oih.opt.jacobianOut.full
                obj.JAll(idxRmv) = [];
            else
                error('jacobianOut must be full or basic!');
            end
            obj.lAll(:,idxRmv) = [];
            if obj.nL>1
                obj.lDir(:,idxRmv) = [];
            end
            obj.signDetJRedAll(:,idxRmv) = [];
            obj.varAll(:,idxRmv) = [];
            if obj.oih.optIsSet.bifAdditionalTestfunction
                obj.bifTestValue(:,idxRmv) = [];
            end
            if obj.oih.stepsizeOptions.iterations
                obj.iterations(:,idxRmv) = [];
            end
            if obj.oih.optIsSet.pathInfoFunction
                obj.pathInfoValue(:,idxRmv) = [];
            end
            if obj.oih.stepsizeOptions.predictor
                obj.xPredictorAll(:,idxRmv) = [];
            end
            if obj.oih.stepsizeOptions.rateOfContraction
                obj.rateOfContraction(:,idxRmv) = [];
            end
            if obj.oih.stepsizeOptions.speedOfContinuation
                obj.speedOfContinuation(:,idxRmv) = [];
            end
        end

        function resetOutput(obj)
            obj.outputFormat = 'singleValueForL';
        end

        function setJacobianByName(obj,J,jacName)
            %% arguments
            arguments
                obj (1,1) continuation.Path
                J (:,:) double
                jacName (1,:) char {mustBeMember(jacName,{'last','previous','initial','plus'})}
            end
            %% get jacobian
            if obj.oih.opt.jacobianOut.basic
                switch jacName
                    case 'initial'
                        obj.JAll{1} = J;
                    case 'last'
                        obj.JAll{2} = J;
                    case 'plus'
                        obj.plusStruct.J = J;
                    case 'previous'
                        obj.JAll{3} = J;
                end
            elseif obj.oih.opt.jacobianOut.full
                switch jacName
                    case 'initial'
                        obj.JAll{1} = J;
                    case 'last'
                        obj.JAll{end} = J;
                    case 'plus'
                        obj.plusStruct.J = J;
                    case 'previous'
                        obj.JAll{end-1} = J;
                end
            else
                error('jacobianOut must be full or basic!');
            end
        end

        function suspend(obj)
            obj.clearPlusStruct();
        end

        function toggleStepback(obj)
            %% arguments
            arguments
                obj (1,1) continuation.Path
            end
            %% check nL
            if obj.nL>1
                obj.stepBackStatus = false;
                obj.oih.do.stepback = false;
                aux.printLine(obj.oih,'------> Stepback must not be called if nL>1.');
            else
                %% toggle plus
                if obj.stepBackStatus
                    %% add plus
                    obj.bifTestValue = [obj.bifTestValue,obj.plusStruct.bifTestValue];
                    if obj.oih.opt.jacobianOut.basic
                        obj.JAll{3} = obj.JAll{2};
                        obj.JAll{2} = obj.plusStruct.J;
                    elseif obj.oih.opt.jacobianOut.full
                        obj.JAll{end+1} = obj.plusStruct.J;
                    else
                        error('jacobianOut must be full or basic!');
                    end
                    obj.lAll = [obj.lAll,obj.plusStruct.l];
                    obj.pathInfoValue = [obj.pathInfoValue,obj.plusStruct.pathInfoValue];
                    obj.signDetJRedAll = [obj.signDetJRedAll,obj.plusStruct.signDetJRed];
                    obj.varAll = [obj.varAll,obj.plusStruct.var];
                    if obj.oih.optIsSet.bifAdditionalTestfunction
                        obj.bifTestValue = [obj.bifTestValue,obj.plusStruct.bifTestValue];
                    end
                    if obj.oih.stepsizeOptions.iterations
                        obj.iterations = [obj.iterations,obj.plusStruct.iterations];
                    end
                    if obj.oih.optIsSet.pathInfoFunction
                        obj.pathInfoValue = [obj.pathInfoValue,obj.plusStruct.pathInfoValue];
                    end
                    if obj.oih.stepsizeOptions.predictor
                        obj.xPredictorAll = [obj.xPredictorAll,obj.plusStruct.xPredictor];
                    end
                    if obj.oih.stepsizeOptions.rateOfContraction
                        obj.rateOfContraction = [obj.rateOfContraction,obj.plusStruct.rateOfContraction];
                    end
                    if obj.oih.stepsizeOptions.speedOfContinuation
                        obj.speedOfContinuation = [obj.speedOfContinuation,obj.plusStruct.speedOfContinuation];
                    end
                    %% clear plus struct
                    obj.clearPlusStruct();
                    %% toggle plus
                    obj.stepBackStatus = false;
                else
                    if obj.nAll<2
                        error("Path must have at least two entries to use 'plus'.");
                    end
                    %% fill plus
                    if obj.oih.opt.jacobianOut.basic
                        obj.plusStruct.J = obj.JAll{2};
                    elseif obj.oih.opt.jacobianOut.full
                        obj.plusStruct.J = obj.JAll{end};
                    else
                        error('jacobianOut must be full or basic!');
                    end
                    obj.plusStruct.l = obj.lAll(:,end);
                    obj.plusStruct.signDetJRed = obj.signDetJRedAll(:,end);
                    obj.plusStruct.var = obj.varAll(:,end);
                    if obj.oih.optIsSet.bifAdditionalTestfunction
                        obj.plusStruct.bifTestValue = obj.bifTestValue(:,end);
                    end
                    if obj.oih.stepsizeOptions.iterations
                        obj.plusStruct.iterations = obj.iterations(:,end);
                    end
                    if obj.oih.optIsSet.pathInfoFunction
                        obj.plusStruct.pathInfoValue = obj.pathInfoValue(:,end);
                    end
                    if obj.oih.stepsizeOptions.predictor
                        obj.plusStruct.xPredictor = obj.xPredictorAll(:,end);
                    end
                    if obj.oih.stepsizeOptions.rateOfContraction
                        obj.plusStruct.rateOfContraction = obj.rateOfContraction(:,end);
                    end
                    if obj.oih.stepsizeOptions.speedOfContinuation
                        obj.plusStruct.speedOfContinuation = obj.speedOfContinuation(:,end);
                    end
                    %% cut all
                    if obj.oih.opt.jacobianOut.basic
                        obj.JAll{2} = obj.JAll{3};
                        obj.JAll{3} = [];
                    elseif obj.oih.opt.jacobianOut.full
                        obj.JAll = obj.JAll(1:(end-1));
                    else
                        error('jacobianOut must be full or basic!');
                    end
                    obj.lAll = obj.lAll(:,1:(end-1));
                    obj.signDetJRedAll = obj.signDetJRedAll(:,1:(end-1));
                    obj.varAll = obj.varAll(:,1:(end-1));
                    if obj.oih.optIsSet.bifAdditionalTestfunction
                        obj.bifTestValue = obj.bifTestValue(:,1:(end-1));
                    end
                    if obj.oih.stepsizeOptions.iterations
                        obj.iterations = obj.iterations(:,1:(end-1));
                    end
                    if obj.oih.optIsSet.pathInfoFunction
                        obj.pathInfoValue = obj.pathInfoValue(:,1:(end-1));
                    end
                    if obj.oih.stepsizeOptions.predictor
                        obj.xPredictorAll = obj.xPredictorAll(:,1:(end-1));
                    end
                    if obj.oih.stepsizeOptions.rateOfContraction
                        obj.rateOfContraction = obj.rateOfContraction(:,1:(end-1));
                    end
                    if obj.oih.stepsizeOptions.speedOfContinuation
                        obj.speedOfContinuation = obj.speedOfContinuation(:,1:(end-1));
                    end
                    %% toggle plus
                    obj.stepBackStatus = true;
                end
            end
        end

        function turn(obj)
            %% arguments
            arguments
                obj (1,1) continuation.Path
            end
            %% deactivate stepback
            obj.stepBackStatus = false;
            obj.clearPlusStruct();
            %% turn path
            idxTurn = obj.nAll:-1:1;
            if obj.oih.opt.jacobianOut.basic
                JLastTemp = obj.JAll{2};
                obj.JAll{2} = obj.JAll{1};
                obj.JAll{3} = [];
                obj.JAll{1} = JLastTemp;
            elseif obj.oih.opt.jacobianOut.full
                obj.JAll = obj.JAll(idxTurn);
            else
                error('jacobianOut must be full or basic!');
            end
            obj.lAll = obj.lAll(:,idxTurn);
            obj.signDetJRedAll = obj.signDetJRedAll(:,idxTurn);
            obj.varAll = obj.varAll(:,idxTurn);
            if obj.oih.optIsSet.bifAdditionalTestfunction
                obj.bifTestValue = obj.bifTestValue(:,idxTurn);
            end
            if obj.oih.stepsizeOptions.iterations
                obj.iterations = obj.iterations(:,idxTurn);
            end
            if obj.oih.optIsSet.pathInfoFunction
                obj.pathInfoValue = obj.pathInfoValue(:,idxTurn);
            end
            if obj.oih.stepsizeOptions.predictor
                obj.xPredictorAll = obj.xPredictorAll(:,idxTurn);
            end
            if obj.oih.stepsizeOptions.rateOfContraction
                obj.rateOfContraction = obj.rateOfContraction(:,idxTurn);
            end
            if obj.oih.stepsizeOptions.speedOfContinuation
                obj.speedOfContinuation = obj.speedOfContinuation(:,idxTurn);
            end
        end

        function vari = varIndex(obj,index)
            arguments
                obj (1,1) continuation.Path
                index (1,:) double {mustBeInteger,mustBeGreaterThanOrEqual(index,1)}
            end
            validateattributes(index,{'double'},'<=',obj.nAll);
            vari = obj.varAll(:,index);
        end

        function xi = xIndex(obj,index)
            arguments
                obj (1,1) continuation.Path
                index (1,:) double {mustBeInteger,mustBeGreaterThanOrEqual(index,1)}
            end
            validateattributes(index,{'double'},'<=',obj.nAll);
            xi = obj.xAll(:,index);
        end
    end

    %% private methods
    methods (Access = private)
        function checkInputOptions(obj,nvaStruct)
            if obj.oih.optIsSet.pathInfoFunction && ~isfield(nvaStruct,'pathInfoValue')
                error('A pathInfoValue must be passed if a pathInfoFunction is set.')
            end
            if obj.oih.optIsSet.bifAdditionalTestfunction && ~isfield(nvaStruct,'bifTestValue')
                error('A bifTestValue must be passed if a bifAdditionalTestfunction is set.')
            end
            if obj.oih.stepsizeOptions.predictor && ~isfield(nvaStruct,'predictor')
                error('predictor must be passed if predictor is set.')
            end
            if obj.oih.stepsizeOptions.speedOfContinuation && ~isfield(nvaStruct,'speedOfContinuation')
                error('speedOfContinuation must be passed if speedOfContinuation is set.')
            end
        end

        function clearPath(obj)
            obj.bifTestValue = [];
            obj.iterations = [];
            obj.JAll = {};
            obj.lAll = [];
            obj.pathInfoValue = [];
            obj.rateOfContraction = [];
            obj.signDetJRedAll = [];
            obj.speedOfContinuation = [];
            obj.varAll = [];
            obj.xPredictorAll = [];
            obj.clearPlusStruct();
        end

        function clearPlusStruct(obj)
            obj.plusStruct.bifTestValue = [];
            obj.plusStruct.iterations = [];
            obj.plusStruct.J = [];
            obj.plusStruct.l = [];
            obj.plusStruct.pathInfoValue = [];
            obj.plusStruct.rateOfContraction = [];
            obj.plusStruct.signDetJRed = [];
            obj.plusStruct.speedOfContinuation = [];
            obj.plusStruct.var = [];
            obj.plusStruct.xPredictor = [];
        end
    end
end