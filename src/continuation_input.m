%% path continuation - continuation_input
%
%   Institute of Dynamics and Vibration Research
%   Leibniz University Hannover
%   08.05.2020 - Alwin Förster
%
function [Opt] = continuation_input(varargin_cell,fun,var0,l_start,l_end)
    %% initialize Opt:
    %
    %  sub-structs:
    %   Sub-structs may only contain true or false values.
    %   Wheter a true value exists can be checkt via 'ison(Opt_sub_struct)'.
    %   If a sub-struct contains multiple true values the first one is valid.
    Opt_homotopy = struct('fix',true,... % st.
                          'newton',false,...
                          'f2',false);
    Opt_sovler = struct('fsolve',true,... % st.
                        'fmincon',false,...
                        'lsqnonlin',false,...
                        'newton',false);
    Opt_arclength = struct('sphere',true,... % st.
                           'linear',false,...
                           'ellipsoid',false,... % coordinate transformed ellipsoid
                           'ellipsoid2',false,... % normal ellipsoid
                           'unique',false); % non-continuation with unique solutions
    Opt_bifurcation = struct('mark',false,...
                             'determine',false,...
                             'trace',false);
    Opt_stepsize = struct('standard',true,...
                          'angle',false,...
                          'curvature',false,...
                          'pid',false);
	%  main-struct:
    Opt = struct('homotopy',Opt_homotopy,...
                 'deflation',true,...
                 'solver',Opt_sovler,...
                 'arclength',Opt_arclength,...
                 'bifurcation',Opt_bifurcation,...
                 'step_size_control',Opt_stepsize,...
                 'step_size_angle',(5/360)*(2*pi),...
                 'step_size_e_max', 0.3,...
                 'step_size_alpha_max', 0.45,...
                 'solver_tol',10^-6,...
                 'display',true,...
                 'max_error_counter',20,...
                 'homotopy_error_counter',10,...
                 'deflation_error_counter',5,...
                 'stepback_error_counter',3,...
                 'jacobian',false,...
                 'l_0',l_start,...
                 'direction',sign(l_end-l_start),...
                 'n_iter_opt',3,...
                 'n_step_max',inf,...
                 'ds_max',inf,...
                 'ds_min',0,...
                 'l_target',l_end,...
                 'alpha_reverse',pi/2,...
                 'unique',false,...
                 'plot',false,...
                 'include_reverse',false,...
                 'predictor_taylor',1,...
                 'predictor_fit',0,...
                 'predictor_adaptive',false,...
                 'plot_vars_index',1:length(var0));
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
                        otherwise
                            error('Unknown parameter %s for option %s',varargin_cell{i+1},varargin_cell{i});
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
                else
                    error('invalid input');
                end
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
    Opt_info = struct('homotopy','#struct',...
                      'deflation','#scalar#logical',...
                      'solver','#struct',...
                      'arclength','#struct',...
                      'bifurcation','#struct',...
                      'step_size_control','#struct',...
                      'step_size_angle','#scalar#double#positive',...
                      'step_size_e_max','#scalar#double#positive',...
                      'step_size_alpha_max','#scalar#double#positive',...
                      'solver_tol','#scalar#double#positive#nonzero',...
                      'display','#scalar#logical',...
                      'max_error_counter','#scalar#integer#positive#nonzero',...
                      'homotopy_error_counter','#scalar#integer#positive#nonzero',...
                      'deflation_error_counter','#scalar#integer#positive#nonzero',...
                      'stepback_error_counter','#scalar#integer#positive#nonzero',...
                      'jacobian','#scalar#logical',...
                      'l_0','#scalar#double',...
                      'direction','#pmone',...
                      'n_iter_opt','#scalar#integer#positive#nonzero',...
                      'n_step_max','#scalar#integer#positive#nonzero',...
                      'ds_max','#scalar#double#positive#nonzero',...
                      'ds_min','#scalar#double#positive#smaller:Opt.ds_max',...
                      'l_target','#scalar#double',...
                      'alpha_reverse','#scalar#double#positive#nonzero',...
                      'unique','#scalar#logical',...
                      'plot','#scalar#logical',...
                      'include_reverse','#scalar#logical',...
                      'predictor_taylor','#scalar#integer#positive#nonzero',...
                      'predictor_fit','#scalar#integer#positive',...
                      'predictor_adaptive','#scalar#logical',...
                      'plot_vars_index','#array#integer#positive#nonzero#unique#max:length(var0)');
    for i=1:numel(Opt_fieldnames)
        if isfield(Opt_info,Opt_fieldnames{i})
            ilower = find(Opt_info.(Opt_fieldnames{i})=='#')+1;
            iupper = [ilower(2:end)-2,length(Opt_info.(Opt_fieldnames{i}))];
            for j=1:length(ilower)
                type = Opt_info.(Opt_fieldnames{i})(ilower(j):iupper(j));
                switch type
                    case 'struct'
                        if ~isa(Opt.(Opt_fieldnames{i}),'struct')
                            error('%s has to be a struct',Opt_fieldnames{i});
                        end
                    case 'logical'
                        if ~isa(Opt.(Opt_fieldnames{i}),'logical')
                            error('%s has to be logical',Opt_fieldnames{i});
                        end
                    case 'double'
                        if ~isa(Opt.(Opt_fieldnames{i}),'double')
                            error('%s has to be double',Opt_fieldnames{i});
                        end
                    case 'integer'
                        if ~prod(floor(Opt.(Opt_fieldnames{i}))==Opt.(Opt_fieldnames{i}))
                            error('%s has to be integer',Opt_fieldnames{i});
                        end
                    case 'positive'
                        if ~prod(Opt.(Opt_fieldnames{i})>=0)
                            error('%s has to be positive',Opt_fieldnames{i});
                        end
                    case 'nonzero'
                        if ~prod(Opt.(Opt_fieldnames{i})~=0)
                            error('%s has to be nonzero',Opt_fieldnames{i});
                        end
                    case 'array'
                        if isempty(Opt.(Opt_fieldnames{i}))
                            error('%s has to be not empty',Opt_fieldnames{i});
                        end
                    case 'scalar'
                        if ~(numel(Opt.(Opt_fieldnames{i}))==1)
                            error('%s has to be scalar',Opt_fieldnames{i});
                        end
                    case 'unique'
                        if ~(numel(Opt.(Opt_fieldnames{i}))==numel(unique(Opt.(Opt_fieldnames{i}))))
                            error('%s has to be unique',Opt_fieldnames{i});
                        end
                    case 'pmone'
                        if ~(Opt.(Opt_fieldnames{i})==+1 || Opt.(Opt_fieldnames{i})==-1)
                            error('%s has to be +-1',Opt_fieldnames{i});
                        end
                    otherwise
                        try
                            if prod(type(1:3)=='max')
                                if ~prod(Opt.(Opt_fieldnames{i})<=eval(type(5:end)))
                                    error('%s has to be smaller equal %.2e',Opt_fieldnames{i},eval(type(5:end)));
                                end
                            elseif prod(type(1:7)=='smaller')
                                if ~prod(Opt.(Opt_fieldnames{i})<=eval(type(9:end)))
                                    error('%s has to be smaller %.2e',Opt_fieldnames{i},eval(type(5:end)));
                                end
                            else
                                error('Field %s has unkown type: %s',Opt_fieldnames{i},type);
                            end
                        catch
                            error('Field %s has unkown type: %s',Opt_fieldnames{i},type);
                        end
                end
            end
        else
            error('Field without type: %s',Opt_fieldnames{i});
        end
    end
    %
end