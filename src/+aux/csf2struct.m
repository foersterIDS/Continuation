%% path continuation - aux.csf2struct
%
%   Institute of Dynamics and Vibration Research
%   Leibniz University Hannover
%   19.08.2021 - Alwin Förster
%
%   Unknown variable names or function handles are capt as strings in the
%   specified field. The output variable evalCell contains the names of
%   fields where this problem occurs. If the variable is known the value of
%   the ith evalCell entry can be reset using the eval function:
%   structOut.(evalCell{i}) = eval(structOut.(evalCell{i}));
%
function [structOut,structInfoOut,evalCell] = csf2struct(structName)
    %% set path
    %
    tempPath = mfilename('fullpath');
    structPath = [tempPath(1:(numel(tempPath)-19)),'settings/',structName,'.csf'];
    clear tempPath;
    %
    %% initialize
    %
    structOut = struct();
    structInfoOut = struct();
    evalCell = {};
    evalCounter = 0;
    sturctFileId = fopen(structPath);
    lineDelimiter = '<>';
    splitDelimiter = '<>';
    lineFormat = ['%s',lineDelimiter,'%s',lineDelimiter,'%s',lineDelimiter,'%s'];
    %
    %% fill struct
    %
    lineStr = fscanf(sturctFileId,lineFormat);
    while ~isempty(lineStr)
        lineCell = strsplit(lineStr,splitDelimiter);
        structInfoOut.(lineCell{1}) = lineCell{4};
        switch lineCell{2}
            case 'eval'
                try
                    structOut.(lineCell{1}) = eval(lineCell{3});
                catch
                    structOut.(lineCell{1}) = lineCell{3};
                    evalCounter = evalCounter+1;
                    evalCell{evalCounter} = lineCell{1};
                end
            case 'double'
                structOut.(lineCell{1}) = str2num(lineCell{3});
            case 'string'
                structOut.(lineCell{1}) = lineCell{3};
            case 'struct'
                tempStruct = aux.csf2struct(lineCell{3});
                structOut.(lineCell{1}) = tempStruct;
                clear tempStruct;
            case 'function_handle'
                structOut.(lineCell{1}) = eval(lineCell{3});
            case 'cell'
                structOut.(lineCell{1}) = eval(lineCell{3});
            otherwise
                error('Unknown field type. Struct-file corrupted?');
        end
        lineStr = fscanf(sturctFileId,lineFormat);
    end
    %
    %% end
    %
    fclose(sturctFileId);
    %
end