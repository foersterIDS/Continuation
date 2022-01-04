%% path continuation - continuation_input
%
%   Institute of Dynamics and Vibration Research
%   Leibniz University Hannover
%   08.05.2020 - Alwin Förster
%   05.01.2021 - Tido Kubatschek
%
function [Opt,ds0] = continuation_input(varargin_cell,fun,var0,l_start,l_end,ds0)
    %% check mandatory input:
    %
    if ~isa(var0,'double') || ~isa(l_start,'double') || ~isa(l_end,'double') || ~isa(ds0,'double')
        error('var0, l_start, l_end and ds0 must be double!');
    end
    if l_start==l_end
        error('l_start and l_end must not be equal!');
    end
    ds0 = abs(ds0);
    %
    %% initialize Opt:
    %
    %  sub-structs:
    %   Sub-structs may only contain true or false values.
    %   Wheter a true value exists can be checked via 'ison(Opt_sub_struct)'.
    %   If a sub-struct contains multiple true values the first one is valid.
    [Opt,Opt_info,eval_cell] = csf2struct('Opt');
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
                                name_legacy = clf2struct('name_legacy');
                            end
                            if isfield(name_legacy.(lower(varargin_cell{i})),lower(varargin_cell{i+1}))
                                varargin_cell{i+1} = name_legacy.(lower(varargin_cell{i})).(lower(varargin_cell{i+1}));
                                i = i-2;
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
                Opt = varargin_cell{i+1};
                for ii=1:numel(Opt_fieldnames)
                    Opt_is_set.(lower(Opt_fieldnames{ii})) = true;
                end
                break;
            else
                %% check name_legacy for Opt-struct:
                if isempty(name_legacy)
                    name_legacy = clf2struct('name_legacy');
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
                                    if (isstruct(Opt.(Opt_fieldnames{i})) && ison(Opt.(Opt_fieldnames{i}))) || (islogical(Opt.(Opt_fieldnames{i})) && Opt.(Opt_fieldnames{i})) || sum(ismember(varargin_cell(1:2:end),Opt_fieldnames{i}))
                                        if isfield(Opt,type(6:end))
                                            fieldname = type(6:end);
                                        else
                                            fieldname = [];
                                        end
                                        try
                                            if ~ison(eval(['Opt.',fieldname]))
                                                if isstruct(Opt.(fieldname))
                                                    Opt_sub_fieldnames = fieldnames(eval(['Opt.',fieldname]));
                                                    eval(['Opt.',fieldname,'.',Opt_sub_fieldnames{1},' = true;']);
                                                    warning('%s has to be on when using %s. %s was set to option %s.',fieldname,Opt_fieldnames{i},fieldname,Opt_sub_fieldnames{1});
                                                elseif islogical(Opt.(fieldname))
                                                    Opt.(fieldname) = true;
                                                    warning('%s has to be on when using %s. %s was set to option ''on''.',fieldname,Opt_fieldnames{i},fieldname);
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
    if ~Opt_is_set.ds_tol
        %% set ds_tol dependent on corrector method:
        if Opt.corrector.sphere
            Opt.ds_tol = [0.99,1.01];
        elseif Opt.corrector.orthogonal
            Opt.ds_tol = [0.99,5];
        elseif Opt.corrector.ellipsoid
            Opt.ds_tol = [0.24,1.01];
        elseif Opt.corrector.ellipsoid2
            Opt.ds_tol = [0.000001,1.01];
        elseif Opt.corrector.unique
            Opt.ds_tol = [0.99,5];
        elseif Opt.corrector.paraboloid
            Opt.ds_tol = [0.09,1.01];
        else
            Opt.ds_tol = [0.5,1.5];
        end
    end
    %
    %% check stepsize options
    %
    if Opt.check_step_size_options
        %
        if ds0 < Opt.ds_min % ds0 must not be lower than ds_min
            error('ds0 cannot be smaller than ds_min. Consider adapting one of them.');
        end
        %
        if ds0 > Opt.ds_max % ds0 must not be larger than ds_max
            error('ds0 cannot be larger than ds_max. Consider adapting one of them.');
        end
        %
        if Opt.ds_min < 10*Opt.solver_tol % ds_min cannot be lower than tolerance of solver
            Opt.ds_min = 10*Opt.solver_tol;
            warning('ds_min has to be at least 10*solver_tol = %.2e. ds_min has been adapted.',Opt.ds_min);
        end
        %
        if ds0 > (sqrt(sum(var0.^2) + (l_end-l_start)^2)/10) % ds_max should not be larger than mag of points
            % get order of magnitude
            n_mag = floor(log10(sqrt(sum(var0.^2) + (l_end-l_start)^2)));
            % adapt ds0
            ds0 = 10^(n_mag-1);
            % order of magnitude must be at least 1 lower
            warning('ds0 is too large! It must not be greater than 1.00e%i. ds0 has been adapted.',n_mag-1);
        end
        %
    end
    %
end