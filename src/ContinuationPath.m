%% path continuation - step_size.StepSizeSingleEvent
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
    end

    %% dependent properties
    properties (Dependent)
        length
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
        
        %% getter:
        function v_all = get.varAll(obj)
            v_all = obj.varAll;
        end

        function l_all = get.lAll(obj)
            l_all = obj.lAll;
        end

        function s_all = get.sAll(obj)
            s_all = obj.sAll;
        end

        function g_all = get.gAll(obj)
            if isempty(obj.gAll) && obj.length>0
                g_all = zeros(0,obj.length);
            else
                g_all = obj.gAll;
            end
        end

        function n = get.length(obj)
            n = numel(obj.lAll);
        end

        function nv = get.nv(obj)
            if obj.length>0
                nv = numel(obj.varAll(:,1));
            else
                nv = NaN;
            end
        end

        function hg = get.hasG(obj)
            if (obj.length>0 && ~isempty(obj.gAll)) || obj.length==0
                hg = true;
            else
                hg = false;
            end
        end

        function xAll = get.xAll(obj)
            xAll = [obj.varAll;obj.lAll;obj.gAll];
        end

        %% setter:
        % ...

        %% misc:
        function addPoint(obj,v,l,g)
            if nargin==3
                if obj.length==0
                    s = 0;
                else
                    dx = obj.xAll(:,end)-[v;l];
                    ds = sqrt(dx'*dx);
                    s = obj.sAll(end)+ds;
                end
                obj.addPointAtS(s,v,l);
            elseif nargin==4
                if obj.length==0
                    s = 0;
                else
                    dx = obj.xAll(:,end)-[v;l;g];
                    ds = sqrt(dx'*dx);
                    s = obj.sAll(end)+ds;
                end
                obj.addPointAtS(s,v,l,'g',g);
            else
                error('Wrong number of input arguments!');
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
            if ~isfield(NameValueArgs,'g') && obj.hasG && obj.length~=0
                error('missing g');
            elseif isfield(NameValueArgs,'g') && ~obj.hasG && obj.length~=0
                error('g cannot be set!');
            end
            if isfield(NameValueArgs,'g')
                g = NameValueArgs.g;
            end
            if obj.length>0
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
            if min(ind)<1 || max(ind)>obj.length
                error('ind must be between 1 and length!');
            end
            CSubPath = ContinuationPath('varAll',obj.varAll(:,ind),'lAll',obj.lAll(:,ind),'gAll',obj.gAll(:,ind),'s0',obj.sAll(ind(1)));
        end

        %% overloaded functions:
        function n = numel(obj)
            n = obj.length();
        end
    end
end