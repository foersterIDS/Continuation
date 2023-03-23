%% path continuation - plot.get_RGB
%
%   Institute of Dynamics and Vibration Research
%   Leibniz University Hannover
%   25.11.2020 - Alwin Förster
%
function [ rgb ] = getRGB( value, maxvalue, minvalue, cm_name )
    value = value(:);
    if nargin<4
        cm_name = 'viridis';
    elseif nargin<3
        error('not enough input arguments');
    elseif nargin>4
        error('to many input arguments');
    end
    cm = colormap(cm_name);
    Pcm = linspace(0,1,numel(cm(:,1)))';
    if maxvalue<minvalue
        temp = minvalue;
        minvalue = maxvalue;
        maxvalue = temp;
    end
    if minvalue==maxvalue
        if value>minvalue
            P = 1;
        elseif value<minvalue
            P = 0;
        else
            P = 0.5;
        end
    else
        P = (value-minvalue)/(maxvalue-minvalue);
        P = max([zeros(size(P)),min([ones(size(P)),P],[],2)],[],2);
    end
    rgb = interp1(Pcm,cm,P);
end