%% path continuation - aux.clf2struct
%
%   Institute of Dynamics and Vibration Research
%   Leibniz University Hannover
%   07.12.2021 - Alwin Förster
%
function [structOut] = clf2struct(structName)
    %% set path
    %
    tempPath = mfilename('fullpath');
    structPath = [tempPath(1:(numel(tempPath)-19)),'settings/',structName,'.clf'];
    clear tempPath;
    %
    %% initialize
    %
    structOut = struct();
    sturctFileId = fopen(structPath);
    lineDelimiter = '>';
    splitDelimiter = '>';
    lineFormat = ['%s',lineDelimiter,'%s',lineDelimiter,'%s'];
    %
    %% fill struct
    %
    lineStr = fscanf(sturctFileId,lineFormat);
    while ~isempty(lineStr)
        lineCell = strsplit(lineStr,splitDelimiter);
        nlc = numel(lineCell);
        if nlc==2
            if isfield(structOut,lineCell{1})
                error('%s struct file corrupted! %s already defined.',structName,lineCell{1});
            else
                structOut.(lineCell{1}) = lineCell{2};
            end
        elseif nlc==3
            if isfield(structOut,lineCell{1})
                if isfield(structOut.(lineCell{1}),lineCell{2})
                    error('%s struct file corrupted! %s.%s already defined.',structName,lineCell{1},lineCell{2});
                else
                    structOut.(lineCell{1}).(lineCell{2}) = lineCell{3};
                end
            else
                structOut.(lineCell{1}).(lineCell{2}) = lineCell{3};
            end
        else
            error('Legacy-file corrupted? (%s)',structName);
        end
        lineStr = fscanf(sturctFileId,lineFormat);
    end
    %
    %% end
    %
    fclose(sturctFileId);
    %
end