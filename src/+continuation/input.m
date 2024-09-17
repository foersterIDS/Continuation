%% path continuation - continuation.input
%
%   Institute of Dynamics and Vibration Research
%   Leibniz University Hannover
%   08.05.2020 - Alwin Förster
%   05.01.2021 - Tido Kubatschek
%
function [oih,ds0,func] = input(vararginCell,fun,var0,lStart,lEnd,ds0)
    %% determine purpose
    %
    purpose = struct('continuation',false,...
                     'homotopy',false,...
                     'parameterTracing',false);
    switch nargin
        case 6
            if abs(nargin(fun))==2
                %% continuation
                % input: vararginCell,fun,var0,lStart,lEnd,ds0
                purpose.continuation = true;
            elseif abs(nargin(fun))==3
                %% parameter tracing
                % input: vararginCell,fun,var0,lStart,lEnd,ds0
                % fun = @(v,l,g) ...
                purpose.parameterTracing = true;
            end
        case 3
            %% homotopy
            % input: vararginCell,fun,var0
            purpose.homotopy = true;
            lStart = 0;
            lEnd = 1;
            ds0 = 0.1;
        otherwise
            error('Illegal use of the function continuation.input(...).');
    end
    %
    %% multi-parameter:
    %
    if numel(lStart)~=1 && ~sum(cellfun(@(x) strcmp(x,'lDirFunction'),vararginCell))
        error('lDirFunction must be set when using multiple parameters.');
    end
    %
    %% check mandatory input:
    %
    if purpose.continuation
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
    [opt,optInfo,evalCell] = aux.csf2struct('Opt');
    for jj=1:numel(evalCell)
        opt.(evalCell{jj}) = eval(opt.(evalCell{jj}));
    end
    OptFieldnames = fieldnames(opt);
    optIsSet = struct();
    for jj=1:numel(OptFieldnames)
        optIsSet.(OptFieldnames{jj}) = false;
        if isstruct(opt.(OptFieldnames{jj}))
            [~,OptSubStructInfo] = aux.csf2struct(['Opt',OptFieldnames{jj}]);
            OptStructInfo.(OptFieldnames{jj}) = OptSubStructInfo;
            OptStructFieldnames.(OptFieldnames{jj}) = fieldnames(OptSubStructInfo);
        end
    end
	%
    %% read vararginCell:
    %
    errMsg = [];
    nameLegacy = [];
    ii = 1;
    while ii<=numel(vararginCell)
        try
            %% adapt to name convention
            OptLowerFieldnames = cellfun(@lower,OptFieldnames,'UniformOutput',false);
            iLower = find(strcmpi(OptLowerFieldnames,vararginCell{ii}));
            if ~isempty(iLower)
                vararginCell{ii} = OptFieldnames{iLower};
            end
            %% set
            if isfield(opt,vararginCell{ii})
                %% set option:
                if ii+1<=numel(vararginCell)
                    if islogical(opt.(vararginCell{ii}))
                        switch vararginCell{ii+1}
                            case 'on'
                                opt.(vararginCell{ii}) = true;
                            case 'off'
                                opt.(vararginCell{ii}) = false;
                            case true
                                opt.(vararginCell{ii}) = true;
                            case false
                                opt.(vararginCell{ii}) = false;
                            otherwise
                                opt.(vararginCell{ii}) = vararginCell{ii+1};
                        end
                    elseif isstruct(opt.(vararginCell{ii}))
                        subOpts = fieldnames(opt.(vararginCell{ii}));
                        %% adapt to name convention
                        OptLowerSubFieldnames = cellfun(@lower,subOpts,'UniformOutput',false);
                        iLower = find(strcmpi(OptLowerSubFieldnames,vararginCell{ii}));
                        if ~isempty(iLower)
                            vararginCell{ii+1} = subOpts{iLower};
                        end
                        %% check
                        if isfield(opt.(vararginCell{ii}),vararginCell{ii+1})
                            for j=1:numel(fieldnames(opt.(vararginCell{ii})))
                                switch vararginCell{ii+1}
                                    case subOpts{j}
                                        opt.(vararginCell{ii}).(subOpts{j}) = true;
                                    otherwise
                                        opt.(vararginCell{ii}).(subOpts{j}) = false;
                                end
                            end
                        elseif strcmpi(vararginCell{ii+1},'on')
                            for j=1:numel(fieldnames(opt.(vararginCell{ii})))
                                switch j
                                    case 1
                                        opt.(vararginCell{ii}).(subOpts{j}) = true;
                                    otherwise
                                        opt.(vararginCell{ii}).(subOpts{j}) = false;
                                end
                            end
                        elseif strcmpi(vararginCell{ii+1},'off')
                            for j=1:numel(fieldnames(opt.(vararginCell{ii})))
                                opt.(vararginCell{ii}).(subOpts{j}) = false;
                            end
                        else
                            %% check nameLegacy for Opt-sub-struct:
                            if isempty(nameLegacy)
                                nameLegacy = aux.clf2struct('nameLegacy');
                            end
                            if isfield(nameLegacy,vararginCell{ii})
                                if isfield(nameLegacy.(vararginCell{ii}),vararginCell{ii+1})
                                    vararginCell{ii+1} = nameLegacy.(vararginCell{ii}).(vararginCell{ii+1});
                                    ii = ii-2;
                                else
                                    errMsg = sprintf('Unknown parameter %s for option %s.',vararginCell{ii+1},vararginCell{ii});
                                    error(errMsg);
                                end
                            else
                                errMsg = sprintf('Unknown parameter %s for option %s.',vararginCell{ii+1},vararginCell{ii});
                                error(errMsg);
                            end
                        end
                    elseif isnumeric(opt.(vararginCell{ii})) && isnumeric(vararginCell{ii+1})
                        opt.(vararginCell{ii}) = vararginCell{ii+1};
                    elseif isa(opt.(vararginCell{ii}),'function_handle') && isa(vararginCell{ii+1},'function_handle')
                        opt.(vararginCell{ii}) = vararginCell{ii+1};
                    elseif iscell(opt.(vararginCell{ii})) && iscell(vararginCell{ii+1})
                        opt.(vararginCell{ii}) = vararginCell{ii+1};
                    else
                        errMsg = sprintf('invalid input');
                        error(errMsg);
                    end
                    optIsSet.(vararginCell{ii}) = true;
                else
                    errMsg = sprintf('Option %s has no value.',vararginCell{ii});
                    error(errMsg);
                end
            elseif strcmpi(vararginCell{ii},'Opt')
                %% set Opt-struct:
                OptTemp = vararginCell{ii+1};
                OptFieldnamesTemp = fieldnames(OptTemp);
                doOptIsSet = true;
                indOptIsSetTemp = find(strcmp(OptFieldnamesTemp,'OptIsSet'));
                if numel(indOptIsSetTemp)>1
                    error('OptIsSet used to often.');
                elseif numel(indOptIsSetTemp)==1
                    optIsSet = OptTemp.(OptFieldnamesTemp{indOptIsSetTemp});
                    OptTemp = rmfield(OptTemp,OptFieldnamesTemp{indOptIsSetTemp});
                    OptFieldnamesTemp(indOptIsSetTemp) = [];
                    doOptIsSet = false;
                end
                usedFields = zeros(numel(OptFieldnamesTemp),1);
                for jj=1:numel(OptFieldnames)
                    containsField = contains(OptFieldnamesTemp,OptFieldnames{jj});
                    if sum(containsField)
                        usedFields(containsField) = 1;
                        opt.(OptFieldnames{jj}) = OptTemp.(OptFieldnames{jj});
                        if doOptIsSet
                            optIsSet.(OptFieldnames{jj}) = true;
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
                    for jj=1:numel(indNotUsed)
                        if jj==1
                            errMsg = [errMsg,OptFieldnamesTemp{indNotUsed(jj)}];
                        elseif jj<numel(indNotUsed)
                            errMsg = [errMsg,', ',OptFieldnamesTemp{indNotUsed(jj)}];
                        else
                            errMsg = [errMsg,' and ',OptFieldnamesTemp{indNotUsed(jj)}];
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
                if isfield(nameLegacy,vararginCell{ii})
                    vararginCell{ii} = nameLegacy.(vararginCell{ii});
                    ii = ii-2;
                else
                    errMsg = sprintf('Unknown option %s.',vararginCell{ii});
                    % list of options:
                    errMsg = [errMsg,' (List of options: ',OptFieldnames{1}];
                    for jj=2:numel(OptFieldnames)
                        errMsg = [errMsg,', ',OptFieldnames{jj}];
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
        ii = ii+2;
    end
    %
    %% read fun:
    %
    if purpose.continuation
        % out:
        if ~(optIsSet.jacobian && ~opt.jacobian)
            try
                if optIsSet.lMult0
                    [R,J] = fun(var0,opt.lMult0);
                else
                    [R,J] = fun(var0,lStart);
                end
                opt.jacobian = true;
            catch
                if optIsSet.lMult0
                    error('A Jacobian must be provided if multiple parameters are tracked.');
                end
                opt.jacobian = false;
            end
        end
        % in:
        if abs(nargin(fun))~=2
            error('%d is an invalid number of input arguments for fun(...).\nfun = fun(v,l)',abs(nargin(fun)))
        else
            func = @(v,l) fun(v,l);
        end
    elseif purpose.parameterTracing
        % out:
        if ~(optIsSet.jacobian && ~opt.jacobian)
            try
                [R,J] = fun(var0,lStart,opt.g0);
                opt.jacobian = true;
            catch
                opt.jacobian = false;
            end
        end
        % in:
        if abs(nargin(fun))~=3
            error('%d is an invalid number of input arguments for fun(...) using parameter tracing.\nfun = fun(v,l,g)',abs(nargin(fun)))
        else
            func = @(v,l) fun(v,l,opt.g0);
        end
        % check settings:
        if ~(opt.bifurcation.parameterTrace || opt.dpa)
            error('fun = @(v,l,g) ... can only be used for DPA.');
        end
    elseif purpose.homotopy
        % out:
        if ~(optIsSet.jacobian && ~opt.jacobian)
            try
                [R,J] = fun(var0);
                opt.jacobian = true;
            catch
                opt.jacobian = false;
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
    [errmsg,opt,optIsSet] = continuation.checkOpt(errmsg,var0,vararginCell,opt,optInfo,OptFieldnames,optIsSet,OptStructInfo,OptStructFieldnames);
    if ~isempty(errmsg)
        errmsg = errmsg(2:end);
        error(errmsg);
    end
    %
    %% set dependent options
    %
    InfoTemp = struct('ds0',ds0,'lEnd',lEnd,'lStart',lStart,'var0',var0);
    opt = aux.updateOpt(opt,optIsSet,InfoTemp);
    %
    %% OptInfoHandle instantiation
    %
    oih = aux.OptInfoHandle(opt,optIsSet);
    %
end