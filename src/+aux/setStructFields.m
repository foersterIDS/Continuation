%% path continuation - aux.setStructFields
%
%   Institute of Dynamics and Vibration Research
%   Leibniz University Hannover
%   10.01.2022 - Tido Kubatschek
%
function [ sStruct ] = setStructFields(sStruct, varargin)
    varLen = nargin-1;
    
    if mod(varLen,2) ~= 0
        error('There must be an even number of fieldnames and values!\n')
    end
    
    for k = 1:2:(varLen-1)
        fieldName = varargin{k};
        
        if ~ischar(fieldName)
            error('Input is not a string!');
        end
        
        fieldValue = varargin{k+1};
        sStruct = setfield(sStruct,fieldName, fieldValue);
    end
    
end

