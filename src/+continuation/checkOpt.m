%% path continuation - continuation.checkOpt
%
%   Institute of Dynamics and Vibration Research
%   Leibniz University Hannover
%   31.03.2022 - Alwin Förster
%
function [errmsg,opt,optIsSet] = checkOpt(errmsg,var0,vararginCell,opt,optInfo,optFieldnames,optIsSet,optStructInfo,optStructFieldnames)
    if nargin>8
        isSubLevel = false;
        optGlobal = opt;
    else
        isSubLevel = true;
        optGlobal = optStructInfo;
    end
    for i=1:numel(optFieldnames)
        if isfield(optInfo,optFieldnames{i})
            icases = find(optInfo.(optFieldnames{i})=='|');
            ncases = length(icases)+1;
            icaseslower = [1,icases+1];
            icasesupper = [icases-1,length(optInfo.(optFieldnames{i}))];
            errmsgTemp = '';
            for k=1:ncases
                errmsgTemp = '';
                types = optInfo.(optFieldnames{i})(icaseslower(k):icasesupper(k));
                ilower = find(types=='#')+1;
                iupper = [ilower(2:end)-2,length(types)];
                for j=1:length(ilower)
                    type = types(ilower(j):iupper(j));
                    switch type
                        case 'array'
                            if isempty(opt.(optFieldnames{i}))
                                errmsgTemp = [errmsgTemp,sprintf('\n%s has to be not empty',optFieldnames{i})];
                            end
                        case 'double'
                            if ~isa(opt.(optFieldnames{i}),'double')
                                errmsgTemp = [errmsgTemp,sprintf('\n%s has to be double',optFieldnames{i})];
                            end
                        case 'false'
                            if isa(opt.(optFieldnames{i}),'logical') && opt.(optFieldnames{i})
                                errmsgTemp = [errmsgTemp,sprintf('\n%s has to be false',optFieldnames{i})];
                            end
                        case 'function_handle'
                            if ~isa(opt.(optFieldnames{i}),'function_handle')
                                errmsgTemp = [errmsgTemp,sprintf('\n%s has to be a function_handle',optFieldnames{i})];
                            end
                        case 'increasing'
                            if ~prod(diff(opt.(optFieldnames{i}))>0)
                                errmsgTemp = [errmsgTemp,sprintf('\n%s must have increasing values',optFieldnames{i})];
                            end
                        case 'integer'
                            if ~prod(floor(opt.(optFieldnames{i}))==opt.(optFieldnames{i}))
                                errmsgTemp = [errmsgTemp,sprintf('\n%s has to be integer',optFieldnames{i})];
                            end
                        case 'isnan'
                            if ~prod(isnan(opt.(optFieldnames{i})))
                                errmsgTemp = [errmsgTemp,sprintf('\n%s has to be NaN',optFieldnames{i})];
                            end
                        case 'logical'
                            if ~isa(opt.(optFieldnames{i}),'logical')
                                errmsgTemp = [errmsgTemp,sprintf('\n%s has to be logical',optFieldnames{i})];
                            end
                        case 'nonzero'
                            if ~prod(opt.(optFieldnames{i})~=0)
                                errmsgTemp = [errmsgTemp,sprintf('\n%s has to be nonzero',optFieldnames{i})];
                            end
                        case 'pmone'
                            if ~(prod(opt.(optFieldnames{i}))==+1 || prod(opt.(optFieldnames{i}))==-1)
                                errmsgTemp = [errmsgTemp,sprintf('\n%s has to be +-1',optFieldnames{i})];
                            end
                        case 'positive'
                            if ~prod(opt.(optFieldnames{i})>=0)
                                errmsgTemp = [errmsgTemp,sprintf('\n%s has to be positive',optFieldnames{i})];
                            end
                        case 'scalar'
                            if ~(numel(opt.(optFieldnames{i}))==1)
                                errmsgTemp = [errmsgTemp,sprintf('\n%s has to be scalar',optFieldnames{i})];
                            end
                        case 'struct'
                            if ~isa(opt.(optFieldnames{i}),'struct')
                                errmsgTemp = [errmsgTemp,sprintf('\n%s has to be a struct',optFieldnames{i})];
                            else
                                [errmsgTemp,opt,optIsSet] = continuation.checkOpt(errmsgTemp,var0,vararginCell,opt.(optFieldnames{i}),optStructInfo.(optFieldnames{i}),optStructFieldnames.(optFieldnames{i}),optIsSet,optGlobal);
                                if isSubLevel
                                    opt = optGlobal;
                                else
                                    optGlobal = opt;
                                end
                            end
                        case 'cell'
                            if ~iscell(opt.(optFieldnames{i}))
                                errmsgTemp = [errmsgTemp,sprintf('\n%s has to be a cell',optFieldnames{i})];
                            end
                        case 'true'
                            if isa(opt.(optFieldnames{i}),'logical') && ~opt.(optFieldnames{i})
                                errmsgTemp = [errmsgTemp,sprintf('\n%s has to be true',optFieldnames{i})];
                            end
                        case 'unique'
                            if ~(numel(opt.(optFieldnames{i}))==numel(unique(opt.(optFieldnames{i}))))
                                errmsgTemp = [errmsgTemp,sprintf('\n%s has to be unique',optFieldnames{i})];
                            end
                        otherwise
                            try
                                if prod(type(1:3)=='max')
                                    if ~prod(opt.(optFieldnames{i})<=eval(type(5:end)))
                                        errmsgTemp = [errmsgTemp,sprintf('\n%s has to be smaller equal %.2e',optFieldnames{i},eval(type(5:end)))];
                                    end
                                elseif prod(type(1:3)=='neq')
                                    if ~prod(opt.(optFieldnames{i})~=eval(type(5:end)))
                                        errmsgTemp = [errmsgTemp,sprintf('\n%s has to be not equal %s',optFieldnames{i},type(5:end))];
                                    end
                                elseif prod(type(1:4)=='ison')
                                    if (isstruct(opt.(optFieldnames{i})) && aux.ison(opt.(optFieldnames{i}))) || (islogical(opt.(optFieldnames{i})) && opt.(optFieldnames{i})) || sum(ismember(vararginCell(1:2:end),optFieldnames{i}))
                                        if isfield(opt,type(6:end))
                                            fieldname = type(6:end);
                                        else
                                            try
                                                splitType = strsplit(type(6:end),'.');
                                                if isfield(optGlobal,splitType{1}) && isfield(optGlobal.(splitType{1}),splitType{2})
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
                                            if ~aux.ison(eval(['optGlobal.',fieldname]))
                                                if isstruct(optGlobal.(fieldname))
                                                    optSubFieldnames = fieldnames(eval(['opt.',fieldname]));
                                                    eval(['opt.',fieldname,'.',optSubFieldnames{1},' = true;']);
                                                    aux.printLine(opt,'--> %s has to be on when using %s. %s was set to option %s.\n',fieldname,optFieldnames{i},fieldname,optSubFieldnames{1});
                                                elseif islogical(opt.(fieldname))
                                                    opt.(fieldname) = true;
                                                    aux.printLine(opt,'--> %s has to be on when using %s. %s was set to option ''on''.\n',fieldname,optFieldnames{i},fieldname);
                                                end
                                            elseif exist('subfieldname')
                                                if ~optGlobal.(fieldname).(subfieldname)
                                                    optGlobal.(fieldname).(subfieldname) = true;
                                                    aux.printLine(opt,'--> %s.%s has to be on when using %s. %s.%s was set to option ''on''.\n',fieldname,subfieldname,optFieldnames{i},fieldname,subfieldname);
                                                    clear subfieldname;
                                                end
                                            end
                                        catch
                                            errmsgTemp = [errmsgTemp,sprintf('\n%s has to be on when using %s',type(6:end),optFieldnames{i})];
                                            clear subfieldname;
                                        end
                                        clear subfieldname;
                                    else
                                        % option not set
                                    end
                                elseif prod(type(1:4)=='norm')
                                    try
                                        if ~(norm(opt.(optFieldnames{i}))==eval(type(6:end)))
                                            opt.(optFieldnames{i}) = opt.(optFieldnames{i})./norm(opt.(optFieldnames{i}));
                                        end
                                    catch
                                        errmsgTemp = [errmsgTemp,sprintf('\n%s has to have norm %.2e but cannot be normed',optFieldnames{i},eval(type(6:end)))];
                                    end
                                elseif prod(type(1:4)=='size')
                                    if ~prod(size(opt.(optFieldnames{i}))==eval(type(6:end)))
                                        errmsgTemp = [errmsgTemp,sprintf('\n%s has to be size [%d,%d]',optFieldnames{i},eval(type(6:end)))];
                                    end
                                elseif prod(type(1:5)=='isoff')
                                    if (isstruct(opt.(optFieldnames{i})) && aux.ison(opt.(optFieldnames{i}))) || (islogical(opt.(optFieldnames{i})) && opt.(optFieldnames{i})) || sum(ismember(vararginCell(1:2:end),optFieldnames{i}))
                                        if (isSubLevel && aux.ison(eval(['optGlobal.',type(7:end)]))) || (~isSubLevel && aux.ison(eval(['opt.',type(7:end)])))
                                            errmsgTemp = [errmsgTemp,sprintf('\n%s has to be off when using %s',type(7:end),optFieldnames{i})];
                                        end
                                    else
                                        % option not set
                                    end
                                elseif prod(type(1:5)=='isset')
                                    if ~optIsSet.(type(7:end)) && (isSubLevel && opt.(optFieldnames{i}))
                                        errmsgTemp = [errmsgTemp,sprintf('\n%s has to be set when using %s',type(7:end),optFieldnames{i})];
                                    end
                                elseif prod(type(1:5)=='seton')
                                    fieldname = type(7:end);
                                    if isSubLevel && opt.(optFieldnames{i})
                                        optGlobal.(fieldname) = true;
                                        optIsSet.(fieldname) = true;
                                    elseif ~isSubLevel
                                        errmsgTemp = [errmsgTemp,sprintf('\nOpt file is corrupted. Incorrect use of seton in option %s.',optFieldnames{i})];
                                    end
                                elseif prod(type(1:6)=='equals')
                                    if ~prod(opt.(optFieldnames{i})==eval(type(8:end)))
                                        errmsgTemp = [errmsgTemp,sprintf('\n%s has to be equal %s',optFieldnames{i},type(8:end))];
                                    end
                                elseif prod(type(1:6)=='larger')
                                    if ~prod(opt.(optFieldnames{i})>eval(type(8:end)))
                                        errmsgTemp = [errmsgTemp,sprintf('\n%s has to be larger %.2e',optFieldnames{i},eval(type(8:end)))];
                                    end
                                elseif prod(type(1:6)=='nargin')
                                    if ~(nargin(opt.(optFieldnames{i}))==eval(type(8:end)))
                                        errmsgTemp = [errmsgTemp,sprintf('\n%s has to be a function_handle with %d inputs',optFieldnames{i},str2num(type(8:end)))];
                                    end
                                elseif prod(type(1:7)=='smaller')
                                    if ~prod(opt.(optFieldnames{i})<=eval(type(9:end)))
                                        errmsgTemp = [errmsgTemp,sprintf('\n%s has to be smaller %.2e',optFieldnames{i},eval(type(9:end)))];
                                    end
                                else
                                    errmsgTemp = [errmsgTemp,sprintf('\nField %s has unkown type: %s',optFieldnames{i},type)];
                                end
                            catch
                                errmsgTemp = [errmsgTemp,sprintf('\nField %s has unkown type: %s',optFieldnames{i},type)];
                            end
                    end
                end
                if isempty(errmsgTemp)
                    break;
                end
            end
            errmsg = [errmsg,errmsgTemp];
        else
            errmsg = [errmsg,sprintf('\nField without type: %s',optFieldnames{i})];
        end
    end
    %% output:
    if isSubLevel
        opt = optGlobal;
    else
        optGlobal = opt;
    end
end