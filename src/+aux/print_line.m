%% path continuation - aux.print_line
%
%   Institute of Dynamics and Vibration Research
%   Leibniz University Hannover
%   19.01.2022 - Alwin Förster
%
function [] = print_line(Opt,str,varargin)
    if Opt.display
        fprintf(str,varargin{:});
    end
end