%% path continuation - aux.csf2struct
%
%   Institute of Dynamics and Vibration Research
%   Leibniz University Hannover
%   19.08.2021 - Alwin Förster
%
%   Unknown variable names or function handles are capt as strings in the
%   specified field. The output variable eval_cell contains the names of
%   fields where this problem occurs. If the variable is known the value of
%   the ith eval_cell entry can be reset using the eval function:
%   struct_out.(eval_cell{i}) = eval(struct_out.(eval_cell{i}));
%
function [struct_out,struct_info_out,eval_cell] = csf2struct(struct_name)
    %% set path
    %
    temp_path = mfilename('fullpath');
    struct_path = [temp_path(1:(numel(temp_path)-19)),'settings\',struct_name,'.csf'];
    clear temp_path;
    %
    %% initialize
    %
    struct_out = struct();
    struct_info_out = struct();
    eval_cell = {};
    eval_counter = 0;
    sturct_file_id = fopen(struct_path);
    line_delimiter = '<>';
    split_delimiter = '<>';
    line_format = ['%s',line_delimiter,'%s',line_delimiter,'%s',line_delimiter,'%s'];
    %
    %% fill struct
    %
    line_str = fscanf(sturct_file_id,line_format);
    while ~isempty(line_str)
        line_cell = strsplit(line_str,split_delimiter);
        struct_info_out.(line_cell{1}) = line_cell{4};
        switch line_cell{2}
            case 'eval'
                try
                    struct_out.(line_cell{1}) = eval(line_cell{3});
                catch
                    struct_out.(line_cell{1}) = line_cell{3};
                    eval_counter = eval_counter+1;
                    eval_cell{eval_counter} = line_cell{1};
                end
            case 'double'
                struct_out.(line_cell{1}) = str2num(line_cell{3});
            case 'string'
                struct_out.(line_cell{1}) = line_cell{3};
            case 'struct'
                temp_struct = aux.csf2struct(line_cell{3});
                struct_out.(line_cell{1}) = temp_struct;
                clear temp_struct;
            case 'function_handle'
                struct_out.(line_cell{1}) = eval(line_cell{3});
            otherwise
                error('Unknown field type. Struct-file corrupted?');
        end
        line_str = fscanf(sturct_file_id,line_format);
    end
    %
    %% end
    %
    fclose(sturct_file_id);
    %
end