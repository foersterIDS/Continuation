%% path continuation - continuation_input
%
%   Institute of Dynamics and Vibration Research
%   Leibniz University Hannover
%   08.05.2020 - Alwin Förster
%   05.01.2020 - Tido Kubatschek
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
    Opt_bifurcation = struct('mark',false,... % st.
                             'determine',false,...
                             'trace',false);
    Opt_corrector = struct('sphere',true,... % st.
                           'linear',false,...
                           'ellipsoid',false,... % coordinate transformed ellipsoid
                           'ellipsoid2',false,... % standard ellipsoid
                           'unique',false); % non-continuation with unique solutions                           
    Opt_diffquot = struct('forward',true,... % st.
                          'central',false);
    Opt_homotopy = struct('fix',true,... % st.
                          'newton',false,...
                          'f2',false);
    Opt_plot = struct('basic',false,... % st.
                      'detail',false,...
                      'semilogx',false,...
                      'semilogy',false,...
                      'loglog',false);              
    Opt_predictor = struct('polynomial',true,... % st.
                           'tangential',false);
	Opt_reverse = struct('angle',true,... % st.
                         'jacobian',false);
	Opt_scaling = struct('dynamicdscale',true,... % st.
                         'staticdscale',false);
    Opt_solver = struct('fsolve',true,... % st.
                        'lsqnonlin',false,...
                        'newton',false);
    Opt_stepsize = struct('iterations',true,... % st.
                          'angle',false,...
                          'curvature',false,...
                          'pid',false);

	%  main-struct:
    Opt = struct('alpha_reverse',pi/2,...
                 'bifurcation',Opt_bifurcation,...
                 'bif_rand_dir', true,...
                 'break_function',@(f,J,v,l,break_fun_out) initial_break_function(f,J,v,l,break_fun_out),...
                 'closed_counter', 3,...
                 'closed_curve_detection',false,...
                 'corrector',Opt_corrector,...
                 'deflation',true,...
                 'deflation_error_counter',5,...
                 'diffquot',Opt_diffquot,...
                 'direction',sign(l_end-l_start)*[zeros(size(var0));1],...
                 'display',true,...
                 'ds_max',inf,...
                 'ds_min',0,...
                 'dscale0',ones(length(var0)+1,1),...
                 'dscale_min', 1e-5,...
                 'homotopy',Opt_homotopy,...
                 'homotopy_error_counter',10,...
                 'include_reverse',false,...
                 'jacobian',false,...
                 'l_0',l_start,...
                 'l_target',l_end,...
                 'live_plot_fig', NaN,...
                 'max_error_counter',20,...
                 'n_bif_search',10,...
                 'n_iter_opt',3,...
                 'n_step_max',inf,...
                 'plot',Opt_plot,...
                 'plot_pause',false,...
                 'plot_var_of_interest', NaN,...
                 'plot_vars_index',1:length(var0),...                
                 'predictor',Opt_predictor,...
                 'predictor_polynomial_adaptive',false,...
                 'predictor_polynomial_fit',0,...
                 'predictor_polynomial_order',1,...
                 'predictor_solver',false,...
                 'reverse',Opt_reverse,...
                 'scaling',Opt_scaling,...
                 'solver',Opt_solver,...
                 'solver_force1it',false,...
                 'solver_tol',10^-6,...
                 'stepback_error_counter',3,...
                 'step_size_alpha_max',0.45,...
                 'step_size_angle',(5/360)*(2*pi),...
                 'step_size_control',Opt_stepsize,...
                 'step_size_e_max',0.3,...
                 'stop_on_bifurcation',false);
	%
    %% read varargin_cell:
    %
    for i=1:2:numel(varargin_cell)-1
        try
            if isfield(Opt,lower(varargin_cell{i}))
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
                        error('Unknown parameter %s for option %s',varargin_cell{i+1},varargin_cell{i});
                    end
                elseif isnumeric(Opt.(lower(varargin_cell{i}))) && isnumeric(varargin_cell{i+1})
                    Opt.(lower(varargin_cell{i})) = varargin_cell{i+1};
                elseif isa(Opt.(lower(varargin_cell{i})),'function_handle') && isa(varargin_cell{i+1},'function_handle')
                    Opt.(lower(varargin_cell{i})) = varargin_cell{i+1};
                else
                    error('invalid input');
                end
            elseif strcmpi(varargin_cell{i},'Opt')
                Opt = varargin_cell{i+1};
                break;
            else
                error('Unknown option %s',varargin_cell{i});
            end
        catch
            error('Check optional input variables');
        end
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
    Opt_fieldnames = fieldnames(Opt);
    Opt_info = struct('alpha_reverse','#scalar#double#positive#nonzero',...
                      'bifurcation','#struct',...
                      'bif_rand_dir','#scalar#logical',...
                      'break_function','#function_handle',...
                      'closed_counter','#scalar#integer#positive#nonzero',...
                      'closed_curve_detection','#scalar#logical',...
                      'corrector','#struct#ison:corrector',...
                      'diffquot','#struct',...
                      'deflation','#scalar#logical',...
                      'deflation_error_counter','#scalar#integer#positive#nonzero',...
                      'direction','#scalar#pmone|#array#double#norm:1#size:[length(var0)+1,1]',...
                      'display','#scalar#logical',...
                      'ds_max','#scalar#double#positive#nonzero',...
                      'ds_min','#scalar#double#positive#smaller:Opt.ds_max',...
                      'dscale0','#array#positive#nonzero#double#size:[length(var0)+1,1]',...
                      'dscale_min','#scalar#positive#nonzero#double|#array#positive#nonzero#double#size:[length(var0)+1,1]',...
                      'homotopy','#struct',...
                      'homotopy_error_counter','#scalar#integer#positive#nonzero',...
                      'include_reverse','#scalar#logical',...
                      'jacobian','#scalar#logical',...
                      'l_0','#scalar#double',...
                      'l_target','#scalar#double',...
                      'live_plot_fig','#scalar#isnan|#scalar#integer#positive#nonzero',...
                      'max_error_counter','#scalar#integer#positive#nonzero',...
                      'n_bif_search','#scalar#integer#positive#nonzero',...
                      'n_iter_opt','#scalar#integer#positive#nonzero',...
                      'n_step_max','#scalar#integer#positive#nonzero',...
                      'plot','#struct',...
                      'plot_pause','#scalar#positive#nonzero#integer|#scalar#logical',...
                      'plot_var_of_interest','#scalar#isnan|#scalar#integer#positive#nonzero#max:length(var0)',...
                      'plot_vars_index','#array#integer#positive#nonzero#unique#max:length(var0)',...
                      'predictor','#struct#ison:predictor',...
                      'predictor_polynomial_adaptive','#scalar#logical',...
                      'predictor_polynomial_fit','#scalar#integer#positive',...
                      'predictor_polynomial_order','#scalar#integer#positive#nonzero',...
                      'predictor_solver','#scalar#logical',...
                      'reverse','#struct',...
                      'scaling','#struct',...
                      'solver','#struct#ison:solver',...
                      'solver_force1it','#scalar#logical',...
                      'solver_tol','#scalar#double#positive#nonzero',...
                      'stepback_error_counter','#scalar#integer#positive#nonzero',...
                      'step_size_alpha_max','#scalar#double#positive#nonzero',...
                      'step_size_angle','#scalar#double#positive',...
                      'step_size_control','#struct#ison:step_size_control',...
                      'step_size_e_max','#scalar#double#positive',...
                      'stop_on_bifurcation','#scalar#logical#ison:bifurcation');
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
                                elseif prod(type(1:4)=='ison')
                                    if (isstruct(Opt.(Opt_fieldnames{i})) || (islogical(Opt.(Opt_fieldnames{i})) && Opt.(Opt_fieldnames{i})))
                                        try
                                            if ~ison(eval(['Opt.',type(6:end)]))
                                                Opt_sub_fieldnames = fieldnames(eval(['Opt.',type(6:end)]));
                                                eval(['Opt.',type(6:end),'.',Opt_sub_fieldnames{1},' = true;']);
                                                warning('%s has to be on. %s was set to option %s.',type(6:end),type(6:end),Opt_sub_fieldnames{1});
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
end