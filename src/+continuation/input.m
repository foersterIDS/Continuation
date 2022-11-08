%% path continuation - continuation.input
%
%   Institute of Dynamics and Vibration Research
%   Leibniz University Hannover
%   08.05.2020 - Alwin Förster
%   05.01.2021 - Tido Kubatschek
%
function [Opt,ds0,Opt_is_set,func] = input(varargin_cell,fun,var0,l_start,l_end,ds0)
    %% determine purpose
    %
    Purpose = struct('continuation',false,...
                     'homotopy',false,...
                     'parameter_tracing',false);
    switch nargin
        case 6
            if abs(nargin(fun))==2
                %% continuation
                % input: varargin_cell,fun,var0,l_start,l_end,ds0
                Purpose.continuation = true;
            elseif abs(nargin(fun))==3
                %% parameter tracing
                % input: varargin_cell,fun,var0,l_start,l_end,ds0
                % fun = @(v,l,g) ...
                Purpose.parameter_tracing = true;
            end
        case 3
            %% homotopy
            % input: varargin_cell,fun,var0
            Purpose.homotopy = true;
            l_start = 0;
            l_end = 1;
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
        if ~isa(var0,'double') || ~isa(l_start,'double') || ~isa(l_end,'double') || ~isa(ds0,'double')
            error('var0, l_start, l_end and ds0 must be double!');
        end
        if l_start==l_end
            error('l_start and l_end must not be equal!');
        end
        ds0 = abs(ds0);
    end
    %
    %% initialize Opt:
    %
    %  sub-structs:
    %   Sub-structs may only contain true or false values.
    %   Wheter a true value exists can be checked via 'aux.ison(Opt_sub_struct)'.
    %   If a sub-struct contains multiple true values the first one is valid.
    [Opt,Opt_info,eval_cell] = aux.csf2struct('Opt');
    for ii=1:numel(eval_cell)
        Opt.(eval_cell{ii}) = eval(Opt.(eval_cell{ii}));
    end
    Opt_fieldnames = fieldnames(Opt);
    Opt_is_set = struct();
    for ii=1:numel(Opt_fieldnames)
        Opt_is_set.(lower(Opt_fieldnames{ii})) = false;
        if isstruct(Opt.(lower(Opt_fieldnames{ii})))
            [~,Opt_sub_struct_info] = aux.csf2struct(['Opt_',Opt_fieldnames{ii}]);
            Opt_struct_info.(lower(Opt_fieldnames{ii})) = Opt_sub_struct_info;
            Opt_struct_fieldnames.(lower(Opt_fieldnames{ii})) = fieldnames(Opt_sub_struct_info);
        end
    end
	%
    %% read varargin_cell:
    %
    err_msg = [];
    name_legacy = [];
    i = 1;
    while i<=numel(varargin_cell)
        try
            if isfield(Opt,lower(varargin_cell{i}))
                %% set option:
                if i+1<=numel(varargin_cell)
                    if islogical(Opt.(lower(varargin_cell{i})))
                        switch lower(varargin_cell{i+1})
                            case 'on'
                                Opt.(lower(varargin_cell{i})) = true;
                            case 'off'
                                Opt.(lower(varargin_cell{i})) = false;
                            case true
                                Opt.(lower(varargin_cell{i})) = true;
                            case false
                                Opt.(lower(varargin_cell{i})) = false;
                            otherwise
                                Opt.(lower(varargin_cell{i})) = varargin_cell{i+1};
                        end
                    elseif isstruct(Opt.(lower(varargin_cell{i})))
                        sub_opts = fieldnames(Opt.(lower(varargin_cell{i})));
                        if isfield(Opt.(lower(varargin_cell{i})),lower(varargin_cell{i+1}))
                            for j=1:numel(fieldnames(Opt.(lower(varargin_cell{i}))))
                                switch lower(varargin_cell{i+1})
                                    case sub_opts{j}
                                        Opt.(lower(varargin_cell{i})).(sub_opts{j}) = true;
                                    otherwise
                                        Opt.(lower(varargin_cell{i})).(sub_opts{j}) = false;
                                end
                            end
                        elseif strcmpi(varargin_cell{i+1},'on')
                            for j=1:numel(fieldnames(Opt.(lower(varargin_cell{i}))))
                                switch j
                                    case 1
                                        Opt.(lower(varargin_cell{i})).(sub_opts{j}) = true;
                                    otherwise
                                        Opt.(lower(varargin_cell{i})).(sub_opts{j}) = false;
                                end
                            end
                        elseif strcmpi(varargin_cell{i+1},'off')
                            for j=1:numel(fieldnames(Opt.(lower(varargin_cell{i}))))
                                Opt.(lower(varargin_cell{i})).(sub_opts{j}) = false;
                            end
                        else
                            %% check name_legacy for Opt-sub-struct:
                            if isempty(name_legacy)
                                name_legacy = aux.clf2struct('name_legacy');
                            end
                            if isfield(name_legacy,lower(varargin_cell{i}))
                                if isfield(name_legacy.(lower(varargin_cell{i})),lower(varargin_cell{i+1}))
                                    varargin_cell{i+1} = name_legacy.(lower(varargin_cell{i})).(lower(varargin_cell{i+1}));
                                    i = i-2;
                                else
                                    err_msg = sprintf('Unknown parameter %s for option %s.',varargin_cell{i+1},varargin_cell{i});
                                    error(err_msg);
                                end
                            else
                                err_msg = sprintf('Unknown parameter %s for option %s.',varargin_cell{i+1},varargin_cell{i});
                                error(err_msg);
                            end
                        end
                    elseif isnumeric(Opt.(lower(varargin_cell{i}))) && isnumeric(varargin_cell{i+1})
                        Opt.(lower(varargin_cell{i})) = varargin_cell{i+1};
                    elseif isa(Opt.(lower(varargin_cell{i})),'function_handle') && isa(varargin_cell{i+1},'function_handle')
                        Opt.(lower(varargin_cell{i})) = varargin_cell{i+1};
                    elseif iscell(Opt.(lower(varargin_cell{i}))) && iscell(varargin_cell{i+1})
                        Opt.(lower(varargin_cell{i})) = varargin_cell{i+1};
                    else
                        err_msg = sprintf('invalid input');
                        error(err_msg);
                    end
                    Opt_is_set.(lower(varargin_cell{i})) = true;
                else
                    err_msg = sprintf('Option %s has no value.',varargin_cell{i});
                    error(err_msg);
                end
            elseif strcmpi(varargin_cell{i},'Opt')
                %% set Opt-struct:
                Opt_temp = varargin_cell{i+1};
                Opt_fieldnames_temp = fieldnames(Opt_temp);
                do_Opt_is_set = true;
                ind_Opt_is_set_temp = find(strcmp(Opt_fieldnames_temp,'Opt_is_set'));
                if numel(ind_Opt_is_set_temp)>1
                    error('Opt_is_set used to often.');
                elseif numel(ind_Opt_is_set_temp)==1
                    Opt_is_set = Opt_temp.(Opt_fieldnames_temp{ind_Opt_is_set_temp});
                    Opt_temp = rmfield(Opt_temp,Opt_fieldnames_temp{ind_Opt_is_set_temp});
                    Opt_fieldnames_temp(ind_Opt_is_set_temp) = [];
                    do_Opt_is_set = false;
                end
                used_fields = zeros(numel(Opt_fieldnames_temp),1);
                for ii=1:numel(Opt_fieldnames)
                    contains_field = contains(Opt_fieldnames_temp,Opt_fieldnames{ii});
                    if sum(contains_field)
                        used_fields(contains_field) = 1;
                        Opt.(lower(Opt_fieldnames{ii})) = Opt_temp.(lower(Opt_fieldnames{ii}));
                        if do_Opt_is_set
                            Opt_is_set.(lower(Opt_fieldnames{ii})) = true;
                        end
                    end
                end
                if ~prod(used_fields)
                    ind_not_used = find(used_fields==0);
                    if numel(ind_not_used)==1
                        err_msg = 'Unknown option ';
                    else
                        err_msg = 'Unknown options ';
                    end
                    for ii=1:numel(ind_not_used)
                        if ii==1
                            err_msg = [err_msg,Opt_fieldnames_temp{ind_not_used(ii)}];
                        elseif ii<numel(ind_not_used)
                            err_msg = [err_msg,', ',Opt_fieldnames_temp{ind_not_used(ii)}];
                        else
                            err_msg = [err_msg,' and ',Opt_fieldnames_temp{ind_not_used(ii)}];
                        end
                    end
                    err_msg = [err_msg,' in user defined Opt-struct.'];
                    error(err_msg);
                end
            else
                %% check name_legacy for Opt-struct:
                if isempty(name_legacy)
                    name_legacy = aux.clf2struct('name_legacy');
                end
                if isfield(name_legacy,lower(varargin_cell{i}))
                    varargin_cell{i} = name_legacy.(lower(varargin_cell{i}));
                    i = i-2;
                else
                    err_msg = sprintf('Unknown option %s.',varargin_cell{i});
                    % list of options:
                    err_msg = [err_msg,' (List of options: ',Opt_fieldnames{1}];
                    for ii=2:numel(Opt_fieldnames)
                        err_msg = [err_msg,', ',Opt_fieldnames{ii}];
                    end
                    error(sprintf(err_msg));
                end
            end
        catch
            if isempty(err_msg)
                error('Unknown error. Check optional input variables.');
            else
                error(err_msg);
            end
        end
        i = i+2;
    end
    %
    %% read fun:
    %
    if Purpose.continuation
        % out:
        try
            [R,J] = fun(var0,l_start);
            Opt.jacobian = true;
        catch
            Opt.jacobian = false;
        end
        % in:
        if abs(nargin(fun))~=2
            error('%d is an invalid number of input arguments for fun(...).\nfun = fun(v,l)',abs(nargin(fun)))
        else
            func = @(v,l) fun(v,l);
        end
    elseif Purpose.parameter_tracing
        % out:
        try
            [R,J] = fun(var0,l_start,Opt.g_0);
            Opt.jacobian = true;
        catch
            Opt.jacobian = false;
        end
        % in:
        if abs(nargin(fun))~=3
            error('%d is an invalid number of input arguments for fun(...) using parameter tracing.\nfun = fun(v,l,g)',abs(nargin(fun)))
        else
            func = @(v,l) fun(v,l,Opt.g_0);
        end
        % check settings:
        if ~(Opt.bifurcation.parameter_trace || Opt.dpa)
            error('fun = @(v,l,g) ... can only be used for DPA.');
        end
    elseif Purpose.homotopy
        % out:
        try
            [R,J] = fun(var0);
            Opt.jacobian = true;
        catch
            Opt.jacobian = false;
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
    [errmsg,Opt,Opt_is_set] = continuation.check_Opt(errmsg,var0,varargin_cell,Opt,Opt_info,Opt_fieldnames,Opt_is_set,Opt_struct_info,Opt_struct_fieldnames);
    if ~isempty(errmsg)
        errmsg = errmsg(2:end);
        error(errmsg);
    end
    %
    %% set dependent options
    %
    Info_temp = struct('ds0',ds0,'l_end',l_end,'l_start',l_start,'var0',var0);
    Opt = aux.update_Opt(Opt,Opt_is_set,Info_temp);
    %
end