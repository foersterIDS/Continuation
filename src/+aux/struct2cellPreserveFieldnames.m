%% path continuation - aux.struct2cellFieldnamePreserving
%
%   Institute of Dynamics and Vibration Research
%   Leibniz University Hannover
%   26.05.2023 - Alwin Förster
%
function c = struct2cellPreserveFieldnames(s)
    arguments
        s (1,1) struct
    end
    fieldNames = fieldnames(s).';
    fieldCell = struct2cell(s).';
    nFields = numel(fieldNames);
    c = cell(1,2*nFields);
    c(1:2:end) = fieldNames;
    c(2:2:end) = fieldCell;
end