%% path continuation - aux.clf2struct
%
%   Institute of Dynamics and Vibration Research
%   Leibniz University Hannover
%   07.12.2021 - Alwin Förster
%
function [struct_out] = clf2struct(struct_name)
    %% set path
    %
    temp_path = mfilename('fullpath');
    struct_path = [temp_path(1:(numel(temp_path)-19)),'settings\',struct_name,'.clf'];
    clear temp_path;
    %
    %% initialize
    %
    struct_out = struct();
    sturct_file_id = fopen(struct_path);
    line_delimiter = '>';
    split_delimiter = '>';
    line_format = ['%s',line_delimiter,'%s',line_delimiter,'%s'];
    %
    %% fill struct
    %
    line_str = fscanf(sturct_file_id,line_format);
    while ~isempty(line_str)
        line_cell = strsplit(line_str,split_delimiter);
        nlc = numel(line_cell);
        if nlc==2
            if isfield(struct_out,line_cell{1})
                error('%s struct file corrupted! %s already defined.',struct_name,line_cell{1});
            else
                struct_out.(line_cell{1}) = line_cell{2};
            end
        elseif nlc==3
            if isfield(struct_out,line_cell{1})
                if isfield(struct_out.(line_cell{1}),line_cell{2})
                    error('%s struct file corrupted! %s.%s already defined.',struct_name,line_cell{1},line_cell{2});
                else
                    struct_out.(line_cell{1}).(line_cell{2}) = line_cell{3};
                end
            else
                struct_out.(line_cell{1}).(line_cell{2}) = line_cell{3};
            end
        else
            error('Legacy-file corrupted? (%s)',struct_name);
        end
        line_str = fscanf(sturct_file_id,line_format);
    end
    %
    %% end
    %
    fclose(sturct_file_id);
    %
end