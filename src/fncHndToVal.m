%% path continuation - fncHndToVal
%
%   Institute of Dynamics and Vibration Research
%   Leibniz University Hannover
%   24.11.2020 - Alwin Förster
%
function [varargout] = fncHndToVal(x,varargin)
    for i=1:nargin-1
        varargout{i} = varargin{i}(x);
    end
end