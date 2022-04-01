%% path continuation - continuation.check_Opt
%
%   Institute of Dynamics and Vibration Research
%   Leibniz University Hannover
%   31.03.2022 - Alwin Förster
%
function [errmsg] = check_Opt(errmsg,var0,varargin_cell,Opt,Opt_info,Opt_fieldnames,Opt_is_set,Opt_struct_info,Opt_struct_fieldnames)
    if nargin>8
        is_sub_level = false;
        Opt_global = Opt;
    else
        is_sub_level = true;
        Opt_global = Opt_struct_info;
    end
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
                        case 'false'
                            if isa(Opt.(Opt_fieldnames{i}),'logical') && Opt.(Opt_fieldnames{i})
                                errmsg_temp = [errmsg_temp,sprintf('\n%s has to be false',Opt_fieldnames{i})];
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
                            else
                                errmsg_temp = continuation.check_Opt(errmsg_temp,var0,varargin_cell,Opt.(Opt_fieldnames{i}),Opt_struct_info.(Opt_fieldnames{i}),Opt_struct_fieldnames.(Opt_fieldnames{i}),Opt_is_set,Opt_global);
                            end
                        case 'true'
                            if isa(Opt.(Opt_fieldnames{i}),'logical') && ~Opt.(Opt_fieldnames{i})
                                errmsg_temp = [errmsg_temp,sprintf('\n%s has to be true',Opt_fieldnames{i})];
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
                                            try
                                                split_type = strsplit(type(6:end),'.');
                                                if isfield(Opt_global,split_type{1}) && isfield(Opt_global.(split_type{1}),split_type{2})
                                                    fieldname = split_type{1};
                                                    subfieldname = split_type{2};
                                                else
                                                    fieldname = [];
                                                end
                                            catch
                                                fieldname = [];
                                                clear subfieldname;
                                            end
                                        end
                                        try
                                            if ~aux.ison(eval(['Opt_global.',fieldname]))
                                                if isstruct(Opt_global.(fieldname))
                                                    Opt_sub_fieldnames = fieldnames(eval(['Opt.',fieldname]));
                                                    eval(['Opt.',fieldname,'.',Opt_sub_fieldnames{1},' = true;']);
                                                    aux.print_line(Opt,'--> %s has to be on when using %s. %s was set to option %s.\n',fieldname,Opt_fieldnames{i},fieldname,Opt_sub_fieldnames{1});
                                                elseif islogical(Opt.(fieldname))
                                                    Opt.(fieldname) = true;
                                                    aux.print_line(Opt,'--> %s has to be on when using %s. %s was set to option ''on''.\n',fieldname,Opt_fieldnames{i},fieldname);
                                                end
                                            elseif exist('subfieldname')
                                                if ~Opt_global.(fieldname).(subfieldname)
                                                    Opt_global.(fieldname).(subfieldname) = true;
                                                    aux.print_line(Opt,'--> %s.%s has to be on when using %s. %s.%s was set to option ''on''.\n',fieldname,subfieldname,Opt_fieldnames{i},fieldname,subfieldname);
                                                    clear subfieldname;
                                                end
                                            end
                                        catch
                                            errmsg_temp = [errmsg_temp,sprintf('\n%s has to be on when using %s',type(6:end),Opt_fieldnames{i})];
                                            clear subfieldname;
                                        end
                                        clear subfieldname;
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
                                elseif prod(type(1:5)=='isoff')
                                    if (isstruct(Opt.(Opt_fieldnames{i})) && aux.ison(Opt.(Opt_fieldnames{i}))) || (islogical(Opt.(Opt_fieldnames{i})) && Opt.(Opt_fieldnames{i})) || sum(ismember(varargin_cell(1:2:end),Opt_fieldnames{i}))
                                        if (is_sub_level && aux.ison(eval(['Opt_global.',type(7:end)]))) || (~is_sub_level && aux.ison(eval(['Opt.',type(7:end)])))
                                            errmsg_temp = [errmsg_temp,sprintf('\n%s has to be off when using %s',type(7:end),Opt_fieldnames{i})];
                                        end
                                    else
                                        % option not set
                                    end
                                elseif prod(type(1:5)=='isset')
                                    if ~Opt_is_set.(type(7:end)) && (is_sub_level && Opt.(Opt_fieldnames{i}))
                                        errmsg_temp = [errmsg_temp,sprintf('\n%s has to be set when using %s',type(7:end),Opt_fieldnames{i})];
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
end