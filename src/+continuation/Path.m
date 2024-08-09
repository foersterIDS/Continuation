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
        Opt
        OptIsSet
        plusStruct = struct('bifTestValue',[],'iterations',[],'J',[],...
                            'pathInfoValue',[],'rateOfContraction',[],...
                            'signDetJ',[],'speedOfContinuation',[],...
                            'var',[],'l',[],'xPredictor',[]);
        StepsizeOptions
    end

    properties (Dependent)
        nAll                    % number of points on path
        sAll                    % arc length
        xAll                    % variables and parameters
    end

    %% public methods
    methods
        %% constructor
        function obj = Path(nVar,nL,Opt,OptIsSet,StepsizeOptions)
            %% arguments
            arguments
                nVar (1,1) double {mustBeInteger,mustBeGreaterThanOrEqual(nVar,1)}
                nL (1,1) double {mustBeInteger,mustBeGreaterThanOrEqual(nL,1)}
                Opt (1,1) struct
                OptIsSet (1,1) struct
                StepsizeOptions (1,1) struct
            end
            %% initialize
            obj.nVar = nVar;
            obj.nL = nL;
            obj.Opt = Opt;
            obj.OptIsSet = OptIsSet;
            obj.StepsizeOptions = StepsizeOptions;
            if obj.Opt.jacobianOut.basic
                obj.JAll = cell(1,3); % initial, last, previous
            elseif obj.Opt.jacobianOut.full
                obj.JAll = {};
            else
                error('jacobianOut must be full or basic!');
            end
            obj.lAll = zeros(nL,0);
            obj.signDetJRedAll = zeros(1,0);
            obj.varAll = zeros(nVar,0);
            obj.resetOutput();
            %% additional functions
            if obj.OptIsSet.pathInfoFunction
                obj.pathInfoValue = zeros(1,0);
            end
            if obj.OptIsSet.bifAdditionalTestfunction
                obj.bifTestValue = zeros(1,0);
            end
            %% step size control values
            if obj.StepsizeOptions.iterations
                obj.iterations = zeros(1,0);
            end
            if obj.StepsizeOptions.predictor
                obj.xPredictorAll = zeros(nVar+nL,0);
            end
            if obj.StepsizeOptions.rateOfContraction
                obj.rateOfContraction = zeros(1,0);
            end
            if obj.StepsizeOptions.speedOfContinuation
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
                        lAll = cumsum([0,sqrt(sum(diff(obj.lAll,1,2).^2,1))]);
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

        function sAll = get.sAll(obj)
            sAll = [0,cumsum(sqrt(sum(diff(obj.xAll,1,2).^2,1)))];
        end

        function xAll = get.xAll(obj)
            xAll = [obj.varAll;obj.lAll];
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
        function addPoint(obj,var,l,J,Solver,pos,NameValueArgs)
            %% arguments
            arguments
                obj (1,1) continuation.Path
                var (:,1) double
                l (:,1) double
                J (:,:) double
                Solver (1,1) struct
                pos (1,1) double {mustBeInteger,mustBeGreaterThan(pos,0)} = obj.nAll+1
                NameValueArgs.bifTestValue (1,1) double
                NameValueArgs.pathInfoValue (1,1) double
                NameValueArgs.predictor (:,1) double
                NameValueArgs.speedOfContinuation (1,1) double
            end
            validateattributes(var,{'numeric'},{'size',[obj.nVar,1]});
            validateattributes(l,{'numeric'},{'size',[obj.nL,1]});
            obj.checkInputOptions(NameValueArgs);
            %% transform parameters
            if obj.nL>1
                if pos~=obj.nAll+1
                    error('With nL>1, entries may only be added at the end of the path.');
                end
                if obj.nAll==0
                    %% initial
                    if numel(l)==obj.nL
                        obj.lRoot = l;
                        if ~obj.OptIsSet.lDirFunction
                            error("'lDirFunction' must be set if more than one parameter is used.");
                        end
                    else
                        error('All parameters must be set at initial point.');
                    end
                else
                    %% additional
                    obj.outputFormat = 'full';
                    lLast = obj.lAll(:,pos-1);
                    obj.resetOutput();
                    if numel(l)==obj.nL
                        %% correct lDir:
                        obj.lDir(:,pos-1) = (l-lLast)/norm(l-lLast);
                    elseif numel(l)==1
                        %% calc. full l
                        lSVAdd = l;
                        lSVEnd = obj.lAll(:,end);
                        dlSV = lSVAdd-lSVEnd;
                        l = lLast+dlSV*obj.lDir(:,end);
                    else
                        error('Number of parameters must be 1 or nL.');
                    end
                end
            end
% #########################################################################
% #########################################################################
% todo: transformation of J??? ############################################
% #########################################################################
% #########################################################################
            %% add/insert
            obj.varAll = [obj.varAll(:,1:(pos-1)),var,obj.varAll(:,pos:end)];
            obj.lAll = [obj.lAll(:,1:(pos-1)),l,obj.lAll(:,pos:end)];
            if obj.nL>1
                lDirTemp = obj.Opt.lDirFunction(var,l,J); % todo #########################################################################################
                lDirTemp = lDirTemp/norm(lDirTemp);
                obj.lDir = [obj.lDir(:,1:(pos-1)),lDirTemp,obj.lDir(:,pos:end)];
            end
            if obj.Opt.jacobianOut.basic
                obj.JAll{3} = obj.JAll{2};
                obj.JAll{2} = J;
            elseif obj.Opt.jacobianOut.full
                obj.JAll = [obj.JAll(:,1:(pos-1)),J,obj.JAll(:,pos:end)];
            else
                error('jacobianOut must be full or basic!');
            end
            obj.signDetJRedAll = [obj.signDetJRedAll(:,1:(pos-1)),sign(det(J(1:obj.nVar,1:obj.nVar))),obj.signDetJRedAll(:,pos:end)];
            if obj.OptIsSet.bifAdditionalTestfunction
                obj.bifTestValue = [obj.bifTestValue(:,1:(pos-1)),NameValueArgs.bifTestValue,obj.bifTestValue(:,pos:end)];
            end
            if obj.StepsizeOptions.iterations
                obj.iterations = [obj.iterations(:,1:(pos-1)),Solver.output.iterations,obj.iterations(:,pos:end)];
            end
            if obj.OptIsSet.pathInfoFunction
                obj.pathInfoValue = [obj.pathInfoValue(:,1:(pos-1)),NameValueArgs.pathInfoValue,obj.pathInfoValue(:,pos:end)];
            end
            if obj.StepsizeOptions.predictor
                obj.xPredictorAll = [obj.xPredictorAll(:,1:(pos-1)),NameValueArgs.predictor,obj.xPredictorAll(:,pos:end)];
            end
            if obj.StepsizeOptions.rateOfContraction
                obj.rateOfContraction = [obj.rateOfContraction(:,1:(pos-1)),Solver.output.rateOfContraction,obj.rateOfContraction(:,pos:end)];
            end
            if obj.StepsizeOptions.speedOfContinuation
                obj.speedOfContinuation = [obj.speedOfContinuation(:,1:(pos-1)),NameValueArgs.speedOfContinuation,obj.speedOfContinuation(:,pos:end)];
            end
        end

        function addPointAtEnd(obj,var,l,J,Solver,NameValueArgs)
            %% arguments
            arguments
                obj (1,1) continuation.Path
                var (:,1) double
                l (:,1) double
                J (:,:) double
                Solver (1,1) struct
                NameValueArgs.bifTestValue (1,1) double
                NameValueArgs.pathInfoValue (1,1) double
                NameValueArgs.predictor (:,1) double
                NameValueArgs.speedOfContinuation (1,1) double
            end
            obj.checkInputOptions(NameValueArgs);
            %% pass to addPoint(...)
            pos = obj.nAll+1;
            nva = aux.struct2cellPreserveFieldnames(NameValueArgs);
            obj.addPoint(var,l,J,Solver,pos,nva{:});
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
            if obj.Opt.jacobianOut.basic
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
            elseif obj.Opt.jacobianOut.full
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
            if obj.Opt.jacobianOut.basic
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
            elseif obj.Opt.jacobianOut.full
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
            if obj.OptIsSet.bifAdditionalTestfunction
                obj.bifTestValue(:,idxRmv) = [];
            end
            if obj.StepsizeOptions.iterations
                obj.iterations(:,idxRmv) = [];
            end
            if obj.OptIsSet.pathInfoFunction
                obj.pathInfoValue(:,idxRmv) = [];
            end
            if obj.StepsizeOptions.predictor
                obj.xPredictorAll(:,idxRmv) = [];
            end
            if obj.StepsizeOptions.rateOfContraction
                obj.rateOfContraction(:,idxRmv) = [];
            end
            if obj.StepsizeOptions.speedOfContinuation
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
            if obj.Opt.jacobianOut.basic
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
            elseif obj.Opt.jacobianOut.full
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
                error('Stepback must not be called if nL>1.');
            end
            %% toggle plus
            if obj.stepBackStatus
                %% add plus
                obj.bifTestValue = [obj.bifTestValue,obj.plusStruct.bifTestValue];
                if obj.Opt.jacobianOut.basic
                    obj.JAll{3} = obj.JAll{2};
                    obj.JAll{2} = obj.plusStruct.J;
                elseif obj.Opt.jacobianOut.full
                    obj.JAll{end+1} = obj.plusStruct.J;
                else
                    error('jacobianOut must be full or basic!');
                end
                obj.lAll = [obj.lAll,obj.plusStruct.l];
                obj.pathInfoValue = [obj.pathInfoValue,obj.plusStruct.pathInfoValue];
                obj.signDetJRedAll = [obj.signDetJRedAll,obj.plusStruct.signDetJRed];
                obj.varAll = [obj.varAll,obj.plusStruct.var];
                if obj.OptIsSet.bifAdditionalTestfunction
                    obj.bifTestValue = [obj.bifTestValue,obj.plusStruct.bifTestValue];
                end
                if obj.StepsizeOptions.iterations
                    obj.iterations = [obj.iterations,obj.plusStruct.iterations];
                end
                if obj.OptIsSet.pathInfoFunction
                    obj.pathInfoValue = [obj.pathInfoValue,obj.plusStruct.pathInfoValue];
                end
                if obj.StepsizeOptions.predictor
                    obj.xPredictorAll = [obj.xPredictorAll,obj.plusStruct.xPredictor];
                end
                if obj.StepsizeOptions.rateOfContraction
                    obj.rateOfContraction = [obj.rateOfContraction,obj.plusStruct.rateOfContraction];
                end
                if obj.StepsizeOptions.speedOfContinuation
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
                if obj.Opt.jacobianOut.basic
                    obj.plusStruct.J = obj.JAll{2};
                elseif obj.Opt.jacobianOut.full
                    obj.plusStruct.J = obj.JAll{end};
                else
                    error('jacobianOut must be full or basic!');
                end
                obj.plusStruct.l = obj.lAll(:,end);
                obj.plusStruct.signDetJRed = obj.signDetJRedAll(:,end);
                obj.plusStruct.var = obj.varAll(:,end);
                if obj.OptIsSet.bifAdditionalTestfunction
                    obj.plusStruct.bifTestValue = obj.bifTestValue(:,end);
                end
                if obj.StepsizeOptions.iterations
                    obj.plusStruct.iterations = obj.iterations(:,end);
                end
                if obj.OptIsSet.pathInfoFunction
                    obj.plusStruct.pathInfoValue = obj.pathInfoValue(:,end);
                end
                if obj.StepsizeOptions.predictor
                    obj.plusStruct.xPredictor = obj.xPredictorAll(:,end);
                end
                if obj.StepsizeOptions.rateOfContraction
                    obj.plusStruct.rateOfContraction = obj.rateOfContraction(:,end);
                end
                if obj.StepsizeOptions.speedOfContinuation
                    obj.plusStruct.speedOfContinuation = obj.speedOfContinuation(:,end);
                end
                %% cut all
                if obj.Opt.jacobianOut.basic
                    obj.JAll{2} = obj.JAll{3};
                    obj.JAll{3} = [];
                elseif obj.Opt.jacobianOut.full
                    obj.JAll = obj.JAll(1:(end-1));
                else
                    error('jacobianOut must be full or basic!');
                end
                obj.lAll = obj.lAll(:,1:(end-1));
                obj.signDetJRedAll = obj.signDetJRedAll(:,1:(end-1));
                obj.varAll = obj.varAll(:,1:(end-1));
                if obj.OptIsSet.bifAdditionalTestfunction
                    obj.bifTestValue = obj.bifTestValue(:,1:(end-1));
                end
                if obj.StepsizeOptions.iterations
                    obj.iterations = obj.iterations(:,1:(end-1));
                end
                if obj.OptIsSet.pathInfoFunction
                    obj.pathInfoValue = obj.pathInfoValue(:,1:(end-1));
                end
                if obj.StepsizeOptions.predictor
                    obj.xPredictorAll = obj.xPredictorAll(:,1:(end-1));
                end
                if obj.StepsizeOptions.rateOfContraction
                    obj.rateOfContraction = obj.rateOfContraction(:,1:(end-1));
                end
                if obj.StepsizeOptions.speedOfContinuation
                    obj.speedOfContinuation = obj.speedOfContinuation(:,1:(end-1));
                end
                %% toggle plus
                obj.stepBackStatus = true;
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
            if obj.Opt.jacobianOut.basic
                JLastTemp = obj.JAll{2};
                obj.JAll{2} = obj.JAll{1};
                obj.JAll{3} = [];
                obj.JAll{1} = JLastTemp;
            elseif obj.Opt.jacobianOut.full
                obj.JAll = obj.JAll(idxTurn);
            else
                error('jacobianOut must be full or basic!');
            end
            obj.lAll = obj.lAll(:,idxTurn);
            obj.signDetJRedAll = obj.signDetJRedAll(:,idxTurn);
            obj.varAll = obj.varAll(:,idxTurn);
            if obj.OptIsSet.bifAdditionalTestfunction
                obj.bifTestValue = obj.bifTestValue(:,idxTurn);
            end
            if obj.StepsizeOptions.iterations
                obj.iterations = obj.iterations(:,idxTurn);
            end
            if obj.OptIsSet.pathInfoFunction
                obj.pathInfoValue = obj.pathInfoValue(:,idxTurn);
            end
            if obj.StepsizeOptions.predictor
                obj.xPredictorAll = obj.xPredictorAll(:,idxTurn);
            end
            if obj.StepsizeOptions.rateOfContraction
                obj.rateOfContraction = obj.rateOfContraction(:,idxTurn);
            end
            if obj.StepsizeOptions.speedOfContinuation
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
            if obj.OptIsSet.pathInfoFunction && ~isfield(nvaStruct,'pathInfoValue')
                error('A pathInfoValue must be passed if a pathInfoFunction is set.')
            end
            if obj.OptIsSet.bifAdditionalTestfunction && ~isfield(nvaStruct,'bifTestValue')
                error('A bifTestValue must be passed if a bifAdditionalTestfunction is set.')
            end
            if obj.StepsizeOptions.predictor && ~isfield(nvaStruct,'predictor')
                error('predictor must be passed if predictor is set.')
            end
            if obj.StepsizeOptions.speedOfContinuation && ~isfield(nvaStruct,'speedOfContinuation')
                error('speedOfContinuation must be passed if speedOfContinuation is set.')
            end
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