clear; close all; clc;
old_dir = pwd;

cd ..\src;
initialfolder = pwd;
listoffolders = dir(initialfolder);
listoffolders = listoffolders(1:end);
dirflags = [listoffolders.isdir];
code_lines = 0;
blank_lines = 0;
comment_lines = 0;
total_lines = 0;
for k = 1 : length(listoffolders)
    if dirflags(k) && ~strcmp(listoffolders(k).name ,'.') && ~strcmp(listoffolders(k).name, '..')
        dirName = listoffolders(k).name;
        subDirs = dir(dirName);
        for k_sub = 3:length(subDirs)
            fileName = subDirs(k_sub).name;
            fullpath = [initialfolder, '\', dirName, '\', fileName];
            [code,comment,blank,total] = count_my_lines(fullpath);
            code_lines = code_lines + code;
            comment_lines = comment_lines + comment;
            blank_lines = blank_lines + blank;
            total_lines = total_lines + total;        
        end
    else
        fileName = listoffolders(k).name;
        if ~strcmp(fileName ,'.') && ~strcmp(fileName, '..')
            [code,comment,blank,total] = count_my_lines(fileName);
            code_lines = code_lines + code;
            comment_lines = comment_lines + comment;
            blank_lines = blank_lines + blank;
            total_lines = total_lines + total;
        end
    end
end

cd(old_dir);
currentday = datetime('now','TimeZone','local','Format','dd-MM-yy');

Lines = [code_lines, comment_lines, blank_lines];

h = figure(); clf;
ax = gca(); 
labels = {'Code: '; 'Comment: '; 'Blank: '};
p = pie(Lines);
newColors = [...
    0, 0.4470, 0.7410;   %blue
    0.4660, 0.6740, 0.1880;   %green
    0.8, 0.8, 0.8];  %gray
ax.Colormap = newColors; 

pText = findobj(p,'Type','text');
percentValues = get(pText,'String'); 
combinedtxt = strcat(labels,percentValues); 
pText(1).String = combinedtxt(1);
pText(2).String = combinedtxt(2);
pText(3).String = combinedtxt(3);
title("Continuation --- Total number of lines: " + string(num2str(total_lines)) + " (" + string(currentday) + ")");
lgd = legend(['Code: ', num2str(code_lines)], ['Comment: ', num2str(comment_lines)], ['Blank: ', num2str(blank_lines)]);
lgd.Location = 'southeastoutside';
h.Color = [1,1,1];
set(h, 'MenuBar', 'none');
set(h, 'ToolBar', 'none');
exportgraphics(h, 'Continuation_Lines.jpg');

function [code,comment,blank,total] = count_my_lines(file)
    fileChar = fileread(file);                                  % read file (char)
    fileStr = strtrim(regexp(fileChar, '\n', 'split'))';        % split lines, remove leading while space
    total = length(fileStr);
    % Remove empty lines
    fileStr(cellfun(@isempty, fileStr)) = [];
    blank = total - length(fileStr);
    % Detect contiguous comments
    startIdx = cumsum(cellfun(@(x)strcmp(x,'%{'), fileStr)); 
    stopIdx = cumsum(cellfun(@(x)strcmp(x,'%}'), fileStr) & startIdx>0);
    contigIdx = (startIdx - stopIdx) > 0; 
    fileStr(contigIdx) = []; 
    % Remove lines that are comments
    fileStr(cellfun(@(x)strcmp(x(1),'%'), fileStr)) = []; 
    % Count number of lines left
    code = length(fileStr);
    comment = total - blank - code;
end