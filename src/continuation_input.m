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
                           'ellipsoid',false);
    Opt_bifurcation = struct('mark',false,...
                             'determine',false,...
                             'trace',false);
	%  main-struct:
    Opt = struct('homotopy',Opt_homotopy,...
                 'deflation',true,...
                 'solver',Opt_sovler,...
                 'arclength',Opt_arclength,...
                 'bifurcation',Opt_bifurcation,...
                 'display',true,...
                 'max_error_counter',10,...
                 'homotopy_error_counter',2,...
                 'jacobian',false,...
                 'l_0',l_start,...
                 'direction',sign(l_end-l_start),...
                 'n_iter_opt',3,...
                 'ds_max',inf);
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
end