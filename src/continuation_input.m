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
    Opt_stepsize = struct('iterations',true,... % st.
                          'angle',false,...
                          'curvature',false,...
                          'pid',false);
    Opt_interpolation = struct('none', true,... %st.
                               'spline', false,...
                               'pchip', false,...
                               'makima', false,...
                               'beziere', false);
	%  main-struct:
    Opt = struct('alpha_reverse',pi/2,...
                 'arclength',Opt_arclength,...
                 'bifurcation',Opt_bifurcation,...
                 'deflation',true,...
                 'deflation_error_counter',5,...
                 'direction',sign(l_end-l_start)*[zeros(size(var0));1],...
                 'display',true,...
                 'ds_max',inf,...
                 'ds_min',0,...
                 'homotopy',Opt_homotopy,...
                 'homotopy_error_counter',10,...
                 'include_reverse',false,...
                 'jacobian',false,...
                 'l_0',l_start,...
                 'l_target',l_end,...
                 'max_error_counter',20,...
                 'n_iter_opt',3,...
                 'n_step_max',inf,...
                 'plot',false,...
                 'plot_vars_index',1:length(var0),...
                 'predictor_adaptive',false,...
                 'predictor_fit',0,...
                 'predictor_taylor',1,...
                 'solver',Opt_sovler,...
                 'solver_tol',10^-6,...
                 'stepback_error_counter',3,...
                 'step_size_alpha_max',0.45,...
                 'step_size_angle',(5/360)*(2*pi),...
                 'step_size_control',Opt_stepsize,...
                 'step_size_e_max',0.3,...
                 'interpolation_method', Opt_interpolation,...
                 'interp_s_inc', 0.01,...
                 'stop_on_bifurcation',false,...
                 'unique',false);
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
    Opt_info = struct('alpha_reverse','#scalar#double#positive#nonzero',...
                      'arclength','#struct',...
                      'bifurcation','#struct',...
                      'deflation','#scalar#logical',...
                      'deflation_error_counter','#scalar#integer#positive#nonzero',...
                      'direction','#scalar#pmone|#array#double#norm:1#size:[length(var0)+1,1]',...
                      'display','#scalar#logical',...
                      'ds_max','#scalar#double#positive#nonzero',...
                      'ds_min','#scalar#double#positive#smaller:Opt.ds_max',...
                      'homotopy','#struct',...
                      'homotopy_error_counter','#scalar#integer#positive#nonzero',...
                      'include_reverse','#scalar#logical',...
                      'jacobian','#scalar#logical',...
                      'l_0','#scalar#double',...
                      'l_target','#scalar#double',...
                      'max_error_counter','#scalar#integer#positive#nonzero',...
                      'n_iter_opt','#scalar#integer#positive#nonzero',...
                      'n_step_max','#scalar#integer#positive#nonzero',...
                      'plot','#scalar#logical',...
                      'plot_vars_index','#array#integer#positive#nonzero#unique#max:length(var0)',...
                      'predictor_adaptive','#scalar#logical',...
                      'predictor_fit','#scalar#integer#positive',...
                      'predictor_taylor','#scalar#integer#positive#nonzero',...
                      'solver','#struct',...
                      'solver_tol','#scalar#double#positive#nonzero',...
                      'stepback_error_counter','#scalar#integer#positive#nonzero',...
                      'step_size_alpha_max','#scalar#double#positive',...
                      'step_size_angle','#scalar#double#positive',...
                      'step_size_control','#struct',...
                      'step_size_e_max','#scalar#double#positive',...
                      'interpolation_method', '#struct',...
                      'interp_s_inc', '#scalar#double#positive#nonzero',...
                      'stop_on_bifurcation','#scalar#logical#ison:bifurcation',...
                      'unique','#scalar#logical');
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
                        case 'integer'
                            if ~prod(floor(Opt.(Opt_fieldnames{i}))==Opt.(Opt_fieldnames{i}))
                                errmsg_temp = [errmsg_temp,sprintf('\n%s has to be integer',Opt_fieldnames{i})];
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
                                    try
                                        if ~ison(eval(['Opt.',type(6:end)]))
                                            Opt_sub_fieldnames = fieldnames(eval(['Opt.',type(6:end)]));
                                            eval(['Opt.',type(6:end),'.',Opt_sub_fieldnames{1},' = true;']);
                                        end
                                    catch
                                        errmsg_temp = [errmsg_temp,sprintf('\n%s has to be on when using %s',type(6:end),Opt_fieldnames{i})];
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
                                elseif prod(type(1:7)=='smaller')
                                    if ~prod(Opt.(Opt_fieldnames{i})<=eval(type(9:end)))
                                        errmsg_temp = [errmsg_temp,sprintf('\n%s has to be smaller %.2e',Opt_fieldnames{i},eval(type(5:end)))];
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