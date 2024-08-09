%% path continuation - aux.printLine
%
%   Institute of Dynamics and Vibration Research
%   Leibniz University Hannover
%   19.01.2022 - Alwin Förster
%
function [] = printLine(oih,str,varargin)
    if oih.opt.display
        fprintf(str,varargin{:});
    end
end