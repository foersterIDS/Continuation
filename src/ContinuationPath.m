%% path continuation - Class: ContinuationPath
%
%   Institute of Dynamics and Vibration Research
%   Leibniz University Hannover
%   08.11.2022 - Alwin Förster
%
classdef ContinuationPath < handle
    %Contains path information
    % v_all : vector of variables (nv x n)
    % l_all : parameter values (1 x n)
    % s_all : arc length (1 x n)
    % (opt)
    % g_all : additional parameter values (1 x n)

    %% public get, private set properties
    properties (GetAccess = public, SetAccess = private)
        varAll
        lAll
        sAll
        gAll
        jacobian
        % TODO:
        % biftest_value
        % bifurcations
        % speed_of_continuation
        % x_predictor
    end

    %% public properties
    properties (Access = public)
        xPredictor
    end

    %% dependent properties
    properties (Dependent)
        ds
        nl
        nv
        hasG
        xAll
    end

    methods
        function obj = ContinuationPath(NameValueArgs)
            arguments
                NameValueArgs.varAll (:,:) double
                NameValueArgs.lAll (1,:) double
                NameValueArgs.gAll (:,:) double
                NameValueArgs.s0 (1,1) double
            end
            if isfield(NameValueArgs,'varAll') && isfield(NameValueArgs,'lAll')
                obj.varAll = NameValueArgs.varAll;
                obj.lAll = NameValueArgs.lAll;
            elseif (isfield(NameValueArgs,'varAll') && ~isfield(NameValueArgs,'lAll')) || (~isfield(NameValueArgs,'varAll') && isfield(NameValueArgs,'lAll'))
                error('If varAll and lAll must both be set.');
            else
                obj.varAll = [];
                obj.lAll = [];
            end
            if isfield(NameValueArgs,'varAll') && isfield(NameValueArgs,'lAll') && isfield(NameValueArgs,'gAll')
                obj.gAll = NameValueArgs.gAll;
            elseif isfield(NameValueArgs,'gAll')
                error('If setting gAll, varAll and lAll must be set.');
            else
                obj.gAll = [];
            end
            if isfield(NameValueArgs,'s0')
                obj.sAll = NameValueArgs.s0;
            else
                obj.sAll = 0;
            end
            obj.updateSAll();
        end



        %% TODO:
        % include_reverse --> addPointAtUknownS
        % get_dscale
        % get_jacobian


        
        %% getter:
        function varAll = get.varAll(obj)
            varAll = obj.varAll;
        end

        function lAll = get.lAll(obj)
            lAll = obj.lAll;
        end

        function sAll = get.sAll(obj)
            sAll = obj.sAll;
        end

        function gAll = get.gAll(obj)
            if isempty(obj.gAll) && obj.nl>0
                gAll = zeros(0,obj.nl);
            else
                gAll = obj.gAll;
            end
        end

        function n = get.nl(obj)
            n = numel(obj.lAll);
        end

        function nv = get.nv(obj)
            if obj.nl>0
                nv = numel(obj.varAll(:,1));
            else
                nv = NaN;
            end
        end

        function hg = get.hasG(obj)
            if (obj.nl>0 && ~isempty(obj.gAll)) || obj.nl==0
                hg = true;
            else
                hg = false;
            end
        end

        function xAll = get.xAll(obj)
            xAll = [obj.varAll;obj.lAll;obj.gAll];
        end

        function ds = get.ds(obj)
            if obj.nl<2
                ds = NaN;
            else
                ds = diff(obj.sAll(end+[-1,0]));
            end
        end

        %% setter:
        % ...

        %% misc:
        function addPoint(obj,v,l,NameValueArgs)
            arguments
                obj (1,1) ContinuationPath
                v (:,1) double
                l (1,1) double
                NameValueArgs.g (1,1) double
                NameValueArgs.Jacobian (:,:) double
            end
            ds = 1;
            if obj.nl==0
                s = 0;
            else
                if ~isfield(NameValueArgs,'g')
                    dx = obj.xAll(:,end)-[v;l];
                else
                    dx = obj.xAll(:,end)-[v;l;NameValueArgs.g];
                end
                ds = sqrt(dx'*dx);
                s = obj.sAll(end)+ds;
            end
            if numel(fieldnames(NameValueArgs))>0
                varArgs = [fieldnames(NameValueArgs).'; struct2cell(NameValueArgs).'];
                obj.addPointAtS(s,v,l,varArgs{:});
            else
                obj.addPointAtS(s,v,l);
            end
            if ds<=0
                error('Point already exists');
            end
        end

        function addPointAtUknownS(obj,v,l,NameValueArgs)
            arguments
                obj (1,1) ContinuationPath
                v (:,1) double
                l (1,1) double
                NameValueArgs.g (1,1) double
            end
            if ~isfield(NameValueArgs,'g') && obj.hasG && obj.nl~=0
                error('missing g');
            elseif isfield(NameValueArgs,'g') && ~obj.hasG && obj.nl~=0
                error('g cannot be set!');
            elseif isfield(NameValueArgs,'g')
                x = [v(:);l;NameValueArgs.g];
            else
                x = [v(:);l];
            end
            if obj.nl>1
                dist = sqrt(sum((obj.xAll-x).^2));
                [d1,i1] = min(dist);
                dist(i1) = inf;
                [d2,i2] = min(dist);
                if abs(i1-i2)==1
                    [ii,jj] = sort([i1,i2]);
                    dd = [d1,d2];
                    dd = dd(jj);
                    s = obj.sAll(ii(1))+dd(1);
                    if isfield(NameValueArgs,'g')
                        obj.addPointAtS(s,v,l,'g',NameValueArgs.g);
                    else
                        obj.addPointAtS(s,v,l);
                    end
                end
            else
                error('Cannot be used for path length <2.');
            end
        end

        function addPointAtS(obj,s,v,l,NameValueArgs)
            arguments
                obj (1,1) ContinuationPath
                s (1,1) double
                v (:,1) double
                l (1,1) double
                NameValueArgs.g (1,1) double
            end
            if ~isfield(NameValueArgs,'g') && obj.hasG && obj.nl~=0
                error('missing g');
            elseif isfield(NameValueArgs,'g') && ~obj.hasG && obj.nl~=0
                error('g cannot be set!');
            end
            if isfield(NameValueArgs,'g')
                g = NameValueArgs.g;
            end
            if obj.nl>0
                [obj.sAll,iSort] = sort([obj.sAll,s]);
                if obj.hasG
                    obj.varAll = [obj.varAll,v(:)];
                    obj.lAll = [obj.lAll,l];
                    obj.gAll = [obj.gAll,g];
                    obj.varAll = obj.varAll(:,iSort);
                    obj.lAll = obj.lAll(iSort);
                    obj.gAll = obj.gAll(iSort);
                else
                    obj.varAll = [obj.varAll,v(:)];
                    obj.lAll = [obj.lAll,l];
                    obj.varAll = obj.varAll(:,iSort);
                    obj.lAll = obj.lAll(iSort);
                end
            else
                obj.varAll = v;
                obj.lAll = l;
                obj.sAll = s;
                if isfield(NameValueArgs,'g')
                    obj.gAll = g;
                else
                    obj.gAll = zeros(0,1);
                end
            end
            obj.updateSAll();
        end

        function updateSAll(obj)
            obj.sAll = [obj.sAll(1),obj.sAll(1)+cumsum(sqrt(sum(diff(obj.xAll,1,2).^2,1)))];
        end

        function CSubPath = copySubPath(obj,ind)
            arguments
                obj (1,1) ContinuationPath
                ind (1,:) double
            end
            if min(ind)<1 || max(ind)>obj.nl
                error('ind must be between 1 and length!');
            end
            CSubPath = ContinuationPath('varAll',obj.varAll(:,ind),'lAll',obj.lAll(:,ind),'gAll',obj.gAll(:,ind),'s0',obj.sAll(ind(1)));
        end

        %% overloaded functions:
        function varargout = numel(obj)
            [varargout{1:nargout}] = obj.nl;
        end

        function varargout = length(obj)
            [varargout{1:nargout}] = obj.nl;
        end
    end
end