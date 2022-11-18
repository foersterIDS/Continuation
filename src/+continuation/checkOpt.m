%% path continuation - continuation.checkOpt
%
%   Institute of Dynamics and Vibration Research
%   Leibniz University Hannover
%   31.03.2022 - Alwin Förster
%
function [errmsg,Opt,OptIsSet] = checkOpt(errmsg,var0,vararginCell,Opt,OptInfo,OptFieldnames,OptIsSet,OptStructInfo,OptStructFieldnames)
    if nargin>8
        isSubLevel = false;
        OptGlobal = Opt;
    else
        isSubLevel = true;
        OptGlobal = OptStructInfo;
    end
    for i=1:numel(OptFieldnames)
        if isfield(OptInfo,OptFieldnames{i})
            icases = find(OptInfo.(OptFieldnames{i})=='|');
            ncases = length(icases)+1;
            icaseslower = [1,icases+1];
            icasesupper = [icases-1,length(OptInfo.(OptFieldnames{i}))];
            errmsgTemp = '';
            for k=1:ncases
                errmsgTemp = '';
                types = OptInfo.(OptFieldnames{i})(icaseslower(k):icasesupper(k));
                ilower = find(types=='#')+1;
                iupper = [ilower(2:end)-2,length(types)];
                for j=1:length(ilower)
                    type = types(ilower(j):iupper(j));
                    switch type
                        case 'array'
                            if isempty(Opt.(OptFieldnames{i}))
                                errmsgTemp = [errmsgTemp,sprintf('\n%s has to be not empty',OptFieldnames{i})];
                            end
                        case 'double'
                            if ~isa(Opt.(OptFieldnames{i}),'double')
                                errmsgTemp = [errmsgTemp,sprintf('\n%s has to be double',OptFieldnames{i})];
                            end
                        case 'false'
                            if isa(Opt.(OptFieldnames{i}),'logical') && Opt.(OptFieldnames{i})
                                errmsgTemp = [errmsgTemp,sprintf('\n%s has to be false',OptFieldnames{i})];
                            end
                        case 'function_handle'
                            if ~isa(Opt.(OptFieldnames{i}),'function_handle')
                                errmsgTemp = [errmsgTemp,sprintf('\n%s has to be a function_handle',OptFieldnames{i})];
                            end
                        case 'increasing'
                            if ~prod(diff(Opt.(OptFieldnames{i}))>0)
                                errmsgTemp = [errmsgTemp,sprintf('\n%s must have increasing values',OptFieldnames{i})];
                            end
                        case 'integer'
                            if ~prod(floor(Opt.(OptFieldnames{i}))==Opt.(OptFieldnames{i}))
                                errmsgTemp = [errmsgTemp,sprintf('\n%s has to be integer',OptFieldnames{i})];
                            end
                        case 'isnan'
                            if ~prod(isnan(Opt.(OptFieldnames{i})))
                                errmsgTemp = [errmsgTemp,sprintf('\n%s has to be NaN',OptFieldnames{i})];
                            end
                        case 'logical'
                            if ~isa(Opt.(OptFieldnames{i}),'logical')
                                errmsgTemp = [errmsgTemp,sprintf('\n%s has to be logical',OptFieldnames{i})];
                            end
                        case 'nonzero'
                            if ~prod(Opt.(OptFieldnames{i})~=0)
                                errmsgTemp = [errmsgTemp,sprintf('\n%s has to be nonzero',OptFieldnames{i})];
                            end
                        case 'pmone'
                            if ~(prod(Opt.(OptFieldnames{i}))==+1 || prod(Opt.(OptFieldnames{i}))==-1)
                                errmsgTemp = [errmsgTemp,sprintf('\n%s has to be +-1',OptFieldnames{i})];
                            end
                        case 'positive'
                            if ~prod(Opt.(OptFieldnames{i})>=0)
                                errmsgTemp = [errmsgTemp,sprintf('\n%s has to be positive',OptFieldnames{i})];
                            end
                        case 'scalar'
                            if ~(numel(Opt.(OptFieldnames{i}))==1)
                                errmsgTemp = [errmsgTemp,sprintf('\n%s has to be scalar',OptFieldnames{i})];
                            end
                        case 'struct'
                            if ~isa(Opt.(OptFieldnames{i}),'struct')
                                errmsgTemp = [errmsgTemp,sprintf('\n%s has to be a struct',OptFieldnames{i})];
                            else
                                [errmsgTemp,Opt,OptIsSet] = continuation.checkOpt(errmsgTemp,var0,vararginCell,Opt.(OptFieldnames{i}),OptStructInfo.(OptFieldnames{i}),OptStructFieldnames.(OptFieldnames{i}),OptIsSet,OptGlobal);
                                if isSubLevel
                                    Opt = OptGlobal;
                                else
                                    OptGlobal = Opt;
                                end
                            end
                        case 'cell'
                            if ~iscell(Opt.(OptFieldnames{i}))
                                errmsgTemp = [errmsgTemp,sprintf('\n%s has to be a cell',OptFieldnames{i})];
                            end
                        case 'true'
                            if isa(Opt.(OptFieldnames{i}),'logical') && ~Opt.(OptFieldnames{i})
                                errmsgTemp = [errmsgTemp,sprintf('\n%s has to be true',OptFieldnames{i})];
                            end
                        case 'unique'
                            if ~(numel(Opt.(OptFieldnames{i}))==numel(unique(Opt.(OptFieldnames{i}))))
                                errmsgTemp = [errmsgTemp,sprintf('\n%s has to be unique',OptFieldnames{i})];
                            end
                        otherwise
                            try
                                if prod(type(1:3)=='max')
                                    if ~prod(Opt.(OptFieldnames{i})<=eval(type(5:end)))
                                        errmsgTemp = [errmsgTemp,sprintf('\n%s has to be smaller equal %.2e',OptFieldnames{i},eval(type(5:end)))];
                                    end
                                elseif prod(type(1:3)=='neq')
                                    if ~prod(Opt.(OptFieldnames{i})~=eval(type(5:end)))
                                        errmsgTemp = [errmsgTemp,sprintf('\n%s has to be not equal %s',OptFieldnames{i},type(5:end))];
                                    end
                                elseif prod(type(1:4)=='ison')
                                    if (isstruct(Opt.(OptFieldnames{i})) && aux.ison(Opt.(OptFieldnames{i}))) || (islogical(Opt.(OptFieldnames{i})) && Opt.(OptFieldnames{i})) || sum(ismember(vararginCell(1:2:end),OptFieldnames{i}))
                                        if isfield(Opt,type(6:end))
                                            fieldname = type(6:end);
                                        else
                                            try
                                                splitType = strsplit(type(6:end),'.');
                                                if isfield(OptGlobal,splitType{1}) && isfield(OptGlobal.(splitType{1}),splitType{2})
                                                    fieldname = splitType{1};
                                                    subfieldname = splitType{2};
                                                else
                                                    fieldname = [];
                                                end
                                            catch
                                                fieldname = [];
                                                clear subfieldname;
                                            end
                                        end
                                        try
                                            if ~aux.ison(eval(['OptGlobal.',fieldname]))
                                                if isstruct(OptGlobal.(fieldname))
                                                    OptSubFieldnames = fieldnames(eval(['Opt.',fieldname]));
                                                    eval(['Opt.',fieldname,'.',OptSubFieldnames{1},' = true;']);
                                                    aux.printLine(Opt,'--> %s has to be on when using %s. %s was set to option %s.\n',fieldname,OptFieldnames{i},fieldname,OptSubFieldnames{1});
                                                elseif islogical(Opt.(fieldname))
                                                    Opt.(fieldname) = true;
                                                    aux.printLine(Opt,'--> %s has to be on when using %s. %s was set to option ''on''.\n',fieldname,OptFieldnames{i},fieldname);
                                                end
                                            elseif exist('subfieldname')
                                                if ~OptGlobal.(fieldname).(subfieldname)
                                                    OptGlobal.(fieldname).(subfieldname) = true;
                                                    aux.printLine(Opt,'--> %s.%s has to be on when using %s. %s.%s was set to option ''on''.\n',fieldname,subfieldname,OptFieldnames{i},fieldname,subfieldname);
                                                    clear subfieldname;
                                                end
                                            end
                                        catch
                                            errmsgTemp = [errmsgTemp,sprintf('\n%s has to be on when using %s',type(6:end),OptFieldnames{i})];
                                            clear subfieldname;
                                        end
                                        clear subfieldname;
                                    else
                                        % option not set
                                    end
                                elseif prod(type(1:4)=='norm')
                                    try
                                        if ~(norm(Opt.(OptFieldnames{i}))==eval(type(6:end)))
                                            Opt.(OptFieldnames{i}) = Opt.(OptFieldnames{i})./norm(Opt.(OptFieldnames{i}));
                                        end
                                    catch
                                        errmsgTemp = [errmsgTemp,sprintf('\n%s has to have norm %.2e but cannot be normed',OptFieldnames{i},eval(type(6:end)))];
                                    end
                                elseif prod(type(1:4)=='size')
                                    if ~prod(size(Opt.(OptFieldnames{i}))==eval(type(6:end)))
                                        errmsgTemp = [errmsgTemp,sprintf('\n%s has to be size [%d,%d]',OptFieldnames{i},eval(type(6:end)))];
                                    end
                                elseif prod(type(1:5)=='isoff')
                                    if (isstruct(Opt.(OptFieldnames{i})) && aux.ison(Opt.(OptFieldnames{i}))) || (islogical(Opt.(OptFieldnames{i})) && Opt.(OptFieldnames{i})) || sum(ismember(vararginCell(1:2:end),OptFieldnames{i}))
                                        if (isSubLevel && aux.ison(eval(['OptGlobal.',type(7:end)]))) || (~isSubLevel && aux.ison(eval(['Opt.',type(7:end)])))
                                            errmsgTemp = [errmsgTemp,sprintf('\n%s has to be off when using %s',type(7:end),OptFieldnames{i})];
                                        end
                                    else
                                        % option not set
                                    end
                                elseif prod(type(1:5)=='isset')
                                    if ~OptIsSet.(type(7:end)) && (isSubLevel && Opt.(OptFieldnames{i}))
                                        errmsgTemp = [errmsgTemp,sprintf('\n%s has to be set when using %s',type(7:end),OptFieldnames{i})];
                                    end
                                elseif prod(type(1:5)=='seton')
                                    fieldname = type(7:end);
                                    if isSubLevel && Opt.(OptFieldnames{i})
                                        OptGlobal.(fieldname) = true;
                                        OptIsSet.(fieldname) = true;
                                    elseif ~isSubLevel
                                        errmsgTemp = [errmsgTemp,sprintf('\nOpt file is corrupted. Incorrect use of seton in option %s.',OptFieldnames{i})];
                                    end
                                elseif prod(type(1:6)=='equals')
                                    if ~prod(Opt.(OptFieldnames{i})==eval(type(8:end)))
                                        errmsgTemp = [errmsgTemp,sprintf('\n%s has to be equal %s',OptFieldnames{i},type(8:end))];
                                    end
                                elseif prod(type(1:6)=='larger')
                                    if ~prod(Opt.(OptFieldnames{i})>eval(type(8:end)))
                                        errmsgTemp = [errmsgTemp,sprintf('\n%s has to be larger %.2e',OptFieldnames{i},eval(type(8:end)))];
                                    end
                                elseif prod(type(1:6)=='nargin')
                                    if ~(nargin(Opt.(OptFieldnames{i}))==eval(type(8:end)))
                                        errmsgTemp = [errmsgTemp,sprintf('\n%s has to be a function_handle with %d inputs',OptFieldnames{i},str2num(type(8:end)))];
                                    end
                                elseif prod(type(1:7)=='smaller')
                                    if ~prod(Opt.(OptFieldnames{i})<=eval(type(9:end)))
                                        errmsgTemp = [errmsgTemp,sprintf('\n%s has to be smaller %.2e',OptFieldnames{i},eval(type(9:end)))];
                                    end
                                else
                                    errmsgTemp = [errmsgTemp,sprintf('\nField %s has unkown type: %s',OptFieldnames{i},type)];
                                end
                            catch
                                errmsgTemp = [errmsgTemp,sprintf('\nField %s has unkown type: %s',OptFieldnames{i},type)];
                            end
                    end
                end
                if isempty(errmsgTemp)
                    break;
                end
            end
            errmsg = [errmsg,errmsgTemp];
        else
            errmsg = [errmsg,sprintf('\nField without type: %s',OptFieldnames{i})];
        end
    end
    %% output:
    if isSubLevel
        Opt = OptGlobal;
    else
        OptGlobal = Opt;
    end
end