%% path continuation - get_RGB
%
%   Institute of Dynamics and Vibration Research
%   Leibniz University Hannover
%   25.11.2020 - Alwin Förster
%
function [ rgb ] = get_RGB( value, maxvalue, minvalue, cm_name )
    if nargin<4
        cm_name = 'parula';
    end
    cm = colormap(cm_name);
    P = (value-minvalue)/(maxvalue-minvalue);
    P = max([0,min([1,P])]);
    Pcm = linspace(0,1,numel(cm(:,1)))';
    rgb = interp1(Pcm,cm,P);
end