%% path continuation - continuation.input
%
%   Institute of Dynamics and Vibration Research
%   Leibniz University Hannover
%   08.05.2020 - Alwin Förster
%   05.01.2021 - Tido Kubatschek
%
function [Opt,ds0,Opt_is_set] = input(varargin_cell,fun,var0,l_start,l_end,ds0)
    %% determine purpose
    %
    Purpose = struct('continuation',false,...
                     'homotopy',false);
    switch nargin
        case 6
            %% continuation
            % input: varargin_cell,fun,var0,l_start,l_end,ds0
            Purpose.continuation = true;
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
                used_fields = zeros(numel(Opt_fieldnames_temp),1);
                for ii=1:numel(Opt_fieldnames)
                    contains_field = contains(Opt_fieldnames_temp,Opt_fieldnames{ii});
                    if sum(contains_field)
                        used_fields(contains_field) = 1;
                        Opt.(lower(Opt_fieldnames{ii})) = Opt_temp.(lower(Opt_fieldnames{ii}));
                        Opt_is_set.(lower(Opt_fieldnames{ii})) = true;
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
                    error(err_msg);
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
        end
    end
    %
    %% check Opt:
    %
	errmsg = '';
    for i=1:numel(Opt_fieldnames)
        if isfield(Opt_info,Opt_fieldnames{i})
            icases = find(Opt_info.(Opt_fieldnames{i})=='|');
            ncases = length(icases)+1;
            icaseslower = [1,icases+1];
            icasesupper = [icases-1,length(Opt_info.(Opt_fieldnames{i}))];
            errmsg_temp = '';
            for k=1:ncases
                errmsg_temp = '';
                types = Opt_info.(Opt_fieldnames{i})(icaseslower(k):icasesupper(k));
                ilower = find(types=='#')+1;
                iupper = [ilower(2:end)-2,length(types)];
                for j=1:length(ilower)
                    type = types(ilower(j):iupper(j));
                    switch type
                        case 'array'
                            if isempty(Opt.(Opt_fieldnames{i}))
                                errmsg_temp = [errmsg_temp,sprintf('\n%s has to be not empty',Opt_fieldnames{i})];
                            end
                        case 'double'
                            if ~isa(Opt.(Opt_fieldnames{i}),'double')
                                errmsg_temp = [errmsg_temp,sprintf('\n%s has to be double',Opt_fieldnames{i})];
                            end
                        case 'function_handle'
                            if ~isa(Opt.(Opt_fieldnames{i}),'function_handle')
                                errmsg_temp = [errmsg_temp,sprintf('\n%s has to be a function_handle',Opt_fieldnames{i})];
                            end
                        case 'increasing'
                            if ~prod(diff(Opt.(Opt_fieldnames{i}))>0)
                                errmsg_temp = [errmsg_temp,sprintf('\n%s must have increasing values',Opt_fieldnames{i})];
                            end
                        case 'integer'
                            if ~prod(floor(Opt.(Opt_fieldnames{i}))==Opt.(Opt_fieldnames{i}))
                                errmsg_temp = [errmsg_temp,sprintf('\n%s has to be integer',Opt_fieldnames{i})];
                            end
                        case 'isnan'
                            if ~prod(isnan(Opt.(Opt_fieldnames{i})))
                                errmsg_temp = [errmsg_temp,sprintf('\n%s has to be NaN',Opt_fieldnames{i})];
                            end
                        case 'logical'
                            if ~isa(Opt.(Opt_fieldnames{i}),'logical')
                                errmsg_temp = [errmsg_temp,sprintf('\n%s has to be logical',Opt_fieldnames{i})];
                            end
                        case 'nonzero'
                            if ~prod(Opt.(Opt_fieldnames{i})~=0)
                                errmsg_temp = [errmsg_temp,sprintf('\n%s has to be nonzero',Opt_fieldnames{i})];
                            end
                        case 'pmone'
                            if ~(prod(Opt.(Opt_fieldnames{i}))==+1 || prod(Opt.(Opt_fieldnames{i}))==-1)
                                errmsg_temp = [errmsg_temp,sprintf('\n%s has to be +-1',Opt_fieldnames{i})];
                            end
                        case 'positive'
                            if ~prod(Opt.(Opt_fieldnames{i})>=0)
                                errmsg_temp = [errmsg_temp,sprintf('\n%s has to be positive',Opt_fieldnames{i})];
                            end
                        case 'scalar'
                            if ~(numel(Opt.(Opt_fieldnames{i}))==1)
                                errmsg_temp = [errmsg_temp,sprintf('\n%s has to be scalar',Opt_fieldnames{i})];
                            end
                        case 'struct'
                            if ~isa(Opt.(Opt_fieldnames{i}),'struct')
                                errmsg_temp = [errmsg_temp,sprintf('\n%s has to be a struct',Opt_fieldnames{i})];
                            end
                        case 'unique'
                            if ~(numel(Opt.(Opt_fieldnames{i}))==numel(unique(Opt.(Opt_fieldnames{i}))))
                                errmsg_temp = [errmsg_temp,sprintf('\n%s has to be unique',Opt_fieldnames{i})];
                            end
                        otherwise
                            try
                                if prod(type(1:3)=='max')
                                    if ~prod(Opt.(Opt_fieldnames{i})<=eval(type(5:end)))
                                        errmsg_temp = [errmsg_temp,sprintf('\n%s has to be smaller equal %.2e',Opt_fieldnames{i},eval(type(5:end)))];
                                    end
                                elseif prod(type(1:3)=='neq')
                                    if ~prod(Opt.(Opt_fieldnames{i})~=eval(type(5:end)))
                                        errmsg_temp = [errmsg_temp,sprintf('\n%s has to be not equal %s',Opt_fieldnames{i},type(5:end))];
                                    end
                                elseif prod(type(1:4)=='ison')
                                    if (isstruct(Opt.(Opt_fieldnames{i})) && aux.ison(Opt.(Opt_fieldnames{i}))) || (islogical(Opt.(Opt_fieldnames{i})) && Opt.(Opt_fieldnames{i})) || sum(ismember(varargin_cell(1:2:end),Opt_fieldnames{i}))
                                        if isfield(Opt,type(6:end))
                                            fieldname = type(6:end);
                                        else
                                            fieldname = [];
                                        end
                                        try
                                            if ~aux.ison(eval(['Opt.',fieldname]))
                                                if isstruct(Opt.(fieldname))
                                                    Opt_sub_fieldnames = fieldnames(eval(['Opt.',fieldname]));
                                                    eval(['Opt.',fieldname,'.',Opt_sub_fieldnames{1},' = true;']);
                                                    aux.print_line(Opt,'--> %s has to be on when using %s. %s was set to option %s.\n',fieldname,Opt_fieldnames{i},fieldname,Opt_sub_fieldnames{1});
                                                elseif islogical(Opt.(fieldname))
                                                    Opt.(fieldname) = true;
                                                    aux.print_line(Opt,'--> %s has to be on when using %s. %s was set to option ''on''.\n',fieldname,Opt_fieldnames{i},fieldname);
                                                end
                                            end
                                        catch
                                            errmsg_temp = [errmsg_temp,sprintf('\n%s has to be on when using %s',type(6:end),Opt_fieldnames{i})];
                                        end
                                    else
                                        % option not set
                                    end
                                elseif prod(type(1:4)=='norm')
                                    try
                                        if ~(norm(Opt.(Opt_fieldnames{i}))==eval(type(6:end)))
                                            Opt.(Opt_fieldnames{i}) = Opt.(Opt_fieldnames{i})./norm(Opt.(Opt_fieldnames{i}));
                                        end
                                    catch
                                        errmsg_temp = [errmsg_temp,sprintf('\n%s has to have norm %.2e but cannot be normed',Opt_fieldnames{i},eval(type(6:end)))];
                                    end
                                elseif prod(type(1:4)=='size')
                                    if ~prod(size(Opt.(Opt_fieldnames{i}))==eval(type(6:end)))
                                        errmsg_temp = [errmsg_temp,sprintf('\n%s has to be size [%d,%d]',Opt_fieldnames{i},eval(type(6:end)))];
                                    end
                                elseif prod(type(1:6)=='equals')
                                    if ~prod(Opt.(Opt_fieldnames{i})==eval(type(8:end)))
                                        errmsg_temp = [errmsg_temp,sprintf('\n%s has to be equal %s',Opt_fieldnames{i},type(8:end))];
                                    end
                                elseif prod(type(1:6)=='larger')
                                    if ~prod(Opt.(Opt_fieldnames{i})>eval(type(8:end)))
                                        errmsg_temp = [errmsg_temp,sprintf('\n%s has to be larger %.2e',Opt_fieldnames{i},eval(type(8:end)))];
                                    end
                                elseif prod(type(1:7)=='smaller')
                                    if ~prod(Opt.(Opt_fieldnames{i})<=eval(type(9:end)))
                                        errmsg_temp = [errmsg_temp,sprintf('\n%s has to be smaller %.2e',Opt_fieldnames{i},eval(type(9:end)))];
                                    end
                                else
                                    errmsg_temp = [errmsg_temp,sprintf('\nField %s has unkown type: %s',Opt_fieldnames{i},type)];
                                end
                            catch
                                errmsg_temp = [errmsg_temp,sprintf('\nField %s has unkown type: %s',Opt_fieldnames{i},type)];
                            end
                    end
                end
                if isempty(errmsg_temp)
                    break;
                end
            end
            errmsg = [errmsg,errmsg_temp];
        else
            errmsg = [errmsg,sprintf('\nField without type: %s',Opt_fieldnames{i})];
        end
    end
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