%% path continuation - continuation.input
%
%   Institute of Dynamics and Vibration Research
%   Leibniz University Hannover
%   08.05.2020 - Alwin Förster
%   05.01.2021 - Tido Kubatschek
%
function [Opt,ds0,OptIsSet,func] = input(vararginCell,fun,var0,lStart,lEnd,ds0)
    %% determine purpose
    %
    Purpose = struct('continuation',false,...
                     'homotopy',false,...
                     'parameterTracing',false);
    switch nargin
        case 6
            if abs(nargin(fun))==2
                %% continuation
                % input: vararginCell,fun,var0,lStart,lEnd,ds0
                Purpose.continuation = true;
            elseif abs(nargin(fun))==3
                %% parameter tracing
                % input: vararginCell,fun,var0,lStart,lEnd,ds0
                % fun = @(v,l,g) ...
                Purpose.parameterTracing = true;
            end
        case 3
            %% homotopy
            % input: vararginCell,fun,var0
            Purpose.homotopy = true;
            lStart = 0;
            lEnd = 1;
            ds0 = 0.1;
        otherwise
            error('Illegal use of the function continuation.input(...).');
    end
    %
    %% check mandatory input:
    %
    if Purpose.continuation
        if isempty(var0)
            error('var0 must not be empty!');
        end
        if ~isa(var0,'double') || ~isa(lStart,'double') || ~isa(lEnd,'double') || ~isa(ds0,'double')
            error('var0, lStart, lEnd and ds0 must be double!');
        end
        if lStart==lEnd
            error('lStart and lEnd must not be equal!');
        end
        ds0 = abs(ds0);
    end
    %
    %% initialize Opt:
    %
    %  sub-structs:
    %   Sub-structs may only contain true or false values.
    %   Wheter a true value exists can be checked via 'aux.ison(OptSubStruct)'.
    %   If a sub-struct contains multiple true values the first one is valid.
    [Opt,OptInfo,evalCell] = aux.csf2struct('Opt');
    for ii=1:numel(evalCell)
        Opt.(evalCell{ii}) = eval(Opt.(evalCell{ii}));
    end
    OptFieldnames = fieldnames(Opt);
    OptIsSet = struct();
    for ii=1:numel(OptFieldnames)
        OptIsSet.(OptFieldnames{ii}) = false;
        if isstruct(Opt.(OptFieldnames{ii}))
            [~,OptSubStructInfo] = aux.csf2struct(['Opt',OptFieldnames{ii}]);
            OptStructInfo.(OptFieldnames{ii}) = OptSubStructInfo;
            OptStructFieldnames.(OptFieldnames{ii}) = fieldnames(OptSubStructInfo);
        end
    end
	%
    %% read vararginCell:
    %
    errMsg = [];
    nameLegacy = [];
    i = 1;
    while i<=numel(vararginCell)
        try
            %% adapt to name convention
            OptLowerFieldnames = cellfun(@lower,OptFieldnames,'UniformOutput',false);
            iLower = find(strcmpi(OptLowerFieldnames,vararginCell{i}));
            if ~isempty(iLower)
                vararginCell{i} = OptFieldnames{iLower};
            end
            %% set
            if isfield(Opt,vararginCell{i})
                %% set option:
                if i+1<=numel(vararginCell)
                    if islogical(Opt.(vararginCell{i}))
                        switch vararginCell{i+1}
                            case 'on'
                                Opt.(vararginCell{i}) = true;
                            case 'off'
                                Opt.(vararginCell{i}) = false;
                            case true
                                Opt.(vararginCell{i}) = true;
                            case false
                                Opt.(vararginCell{i}) = false;
                            otherwise
                                Opt.(vararginCell{i}) = vararginCell{i+1};
                        end
                    elseif isstruct(Opt.(vararginCell{i}))
                        subOpts = fieldnames(Opt.(vararginCell{i}));
                        %% adapt to name convention
                        OptLowerSubFieldnames = cellfun(@lower,subOpts,'UniformOutput',false);
                        iLower = find(strcmpi(OptLowerSubFieldnames,vararginCell{i}));
                        if ~isempty(iLower)
                            vararginCell{i+1} = subOpts{iLower};
                        end
                        %% check
                        if isfield(Opt.(vararginCell{i}),vararginCell{i+1})
                            for j=1:numel(fieldnames(Opt.(vararginCell{i})))
                                switch vararginCell{i+1}
                                    case subOpts{j}
                                        Opt.(vararginCell{i}).(subOpts{j}) = true;
                                    otherwise
                                        Opt.(vararginCell{i}).(subOpts{j}) = false;
                                end
                            end
                        elseif strcmpi(vararginCell{i+1},'on')
                            for j=1:numel(fieldnames(Opt.(vararginCell{i})))
                                switch j
                                    case 1
                                        Opt.(vararginCell{i}).(subOpts{j}) = true;
                                    otherwise
                                        Opt.(vararginCell{i}).(subOpts{j}) = false;
                                end
                            end
                        elseif strcmpi(vararginCell{i+1},'off')
                            for j=1:numel(fieldnames(Opt.(vararginCell{i})))
                                Opt.(vararginCell{i}).(subOpts{j}) = false;
                            end
                        else
                            %% check nameLegacy for Opt-sub-struct:
                            if isempty(nameLegacy)
                                nameLegacy = aux.clf2struct('nameLegacy');
                            end
                            if isfield(nameLegacy,vararginCell{i})
                                if isfield(nameLegacy.(vararginCell{i}),vararginCell{i+1})
                                    vararginCell{i+1} = nameLegacy.(vararginCell{i}).(vararginCell{i+1});
                                    i = i-2;
                                else
                                    errMsg = sprintf('Unknown parameter %s for option %s.',vararginCell{i+1},vararginCell{i});
                                    error(errMsg);
                                end
                            else
                                errMsg = sprintf('Unknown parameter %s for option %s.',vararginCell{i+1},vararginCell{i});
                                error(errMsg);
                            end
                        end
                    elseif isnumeric(Opt.(vararginCell{i})) && isnumeric(vararginCell{i+1})
                        Opt.(vararginCell{i}) = vararginCell{i+1};
                    elseif isa(Opt.(vararginCell{i}),'function_handle') && isa(vararginCell{i+1},'function_handle')
                        Opt.(vararginCell{i}) = vararginCell{i+1};
                    elseif iscell(Opt.(vararginCell{i})) && iscell(vararginCell{i+1})
                        Opt.(vararginCell{i}) = vararginCell{i+1};
                    else
                        errMsg = sprintf('invalid input');
                        error(errMsg);
                    end
                    OptIsSet.(vararginCell{i}) = true;
                else
                    errMsg = sprintf('Option %s has no value.',vararginCell{i});
                    error(errMsg);
                end
            elseif strcmpi(vararginCell{i},'Opt')
                %% set Opt-struct:
                OptTemp = vararginCell{i+1};
                OptFieldnamesTemp = fieldnames(OptTemp);
                doOptIsSet = true;
                indOptIsSetTemp = find(strcmp(OptFieldnamesTemp,'OptIsSet'));
                if numel(indOptIsSetTemp)>1
                    error('OptIsSet used to often.');
                elseif numel(indOptIsSetTemp)==1
                    OptIsSet = OptTemp.(OptFieldnamesTemp{indOptIsSetTemp});
                    OptTemp = rmfield(OptTemp,OptFieldnamesTemp{indOptIsSetTemp});
                    OptFieldnamesTemp(indOptIsSetTemp) = [];
                    doOptIsSet = false;
                end
                usedFields = zeros(numel(OptFieldnamesTemp),1);
                for ii=1:numel(OptFieldnames)
                    containsField = contains(OptFieldnamesTemp,OptFieldnames{ii});
                    if sum(containsField)
                        usedFields(containsField) = 1;
                        Opt.(OptFieldnames{ii}) = OptTemp.(OptFieldnames{ii});
                        if doOptIsSet
                            OptIsSet.(OptFieldnames{ii}) = true;
                        end
                    end
                end
                if ~prod(usedFields)
                    indNotUsed = find(usedFields==0);
                    if numel(indNotUsed)==1
                        errMsg = 'Unknown option ';
                    else
                        errMsg = 'Unknown options ';
                    end
                    for ii=1:numel(indNotUsed)
                        if ii==1
                            errMsg = [errMsg,OptFieldnamesTemp{indNotUsed(ii)}];
                        elseif ii<numel(indNotUsed)
                            errMsg = [errMsg,', ',OptFieldnamesTemp{indNotUsed(ii)}];
                        else
                            errMsg = [errMsg,' and ',OptFieldnamesTemp{indNotUsed(ii)}];
                        end
                    end
                    errMsg = [errMsg,' in user defined Opt-struct.'];
                    error(errMsg);
                end
            else
                %% check nameLegacy for Opt-struct:
                if isempty(nameLegacy)
                    nameLegacy = aux.clf2struct('nameLegacy');
                end
                if isfield(nameLegacy,vararginCell{i})
                    vararginCell{i} = nameLegacy.(vararginCell{i});
                    i = i-2;
                else
                    errMsg = sprintf('Unknown option %s.',vararginCell{i});
                    % list of options:
                    errMsg = [errMsg,' (List of options: ',OptFieldnames{1}];
                    for ii=2:numel(OptFieldnames)
                        errMsg = [errMsg,', ',OptFieldnames{ii}];
                    end
                    error(sprintf(errMsg));
                end
            end
        catch
            if isempty(errMsg)
                error('Unknown error. Check optional input variables.');
            else
                error(errMsg);
            end
        end
        i = i+2;
    end
    %
    %% read fun:
    %
    if Purpose.continuation
        % out:
        if ~(OptIsSet.jacobian && ~Opt.jacobian)
            try
                [R,J] = fun(var0,lStart);
                Opt.jacobian = true;
            catch
                Opt.jacobian = false;
            end
        end
        % in:
        if abs(nargin(fun))~=2
            error('%d is an invalid number of input arguments for fun(...).\nfun = fun(v,l)',abs(nargin(fun)))
        else
            func = @(v,l) fun(v,l);
        end
    elseif Purpose.parameterTracing
        % out:
        if ~(OptIsSet.jacobian && ~Opt.jacobian)
            try
                [R,J] = fun(var0,lStart,Opt.g0);
                Opt.jacobian = true;
            catch
                Opt.jacobian = false;
            end
        end
        % in:
        if abs(nargin(fun))~=3
            error('%d is an invalid number of input arguments for fun(...) using parameter tracing.\nfun = fun(v,l,g)',abs(nargin(fun)))
        else
            func = @(v,l) fun(v,l,Opt.g0);
        end
        % check settings:
        if ~(Opt.bifurcation.parameterTrace || Opt.dpa)
            error('fun = @(v,l,g) ... can only be used for DPA.');
        end
    elseif Purpose.homotopy
        % out:
        if ~(OptIsSet.jacobian && ~Opt.jacobian)
            try
                [R,J] = fun(var0);
                Opt.jacobian = true;
            catch
                Opt.jacobian = false;
            end
        end
        % in:
        if abs(nargin(fun))~=1
            error('%d is an invalid number of input arguments for fun(...).\nfun = fun(v)',abs(nargin(fun)))
        else
            func = @(v) fun(v);
        end
    end
    %
    %% check Opt:
    %
	errmsg = '';
    [errmsg,Opt,OptIsSet] = continuation.checkOpt(errmsg,var0,vararginCell,Opt,OptInfo,OptFieldnames,OptIsSet,OptStructInfo,OptStructFieldnames);
    if ~isempty(errmsg)
        errmsg = errmsg(2:end);
        error(errmsg);
    end
    %
    %% set dependent options
    %
    InfoTemp = struct('ds0',ds0,'lEnd',lEnd,'lStart',lStart,'var0',var0);
    Opt = aux.updateOpt(Opt,OptIsSet,InfoTemp);
    %
end