%% path continuation - plot.getRGB
%
%   Institute of Dynamics and Vibration Research
%   Leibniz University Hannover
%   25.11.2020 - Alwin Förster
%
function [ rgb ] = getRGB( value, maxValue, minValue, cmName )
    value = value(:);
    if nargin==1 && isa(value,'char')
        switch value
            case 'r'
                rgb = [183,39,58]/255;
            case 'b'
                rgb = [65,76,204]/255;
            case 'g'
                rgb = [65,163,52]/255;
            case 'y'
                rgb = [183,142,39]/255;
            case 'e'
                rgb = [181,0,24]/255;
            case 'p'
                rgb = [0,20,204]/255;
            case 'k'
                rgb = [0,0,0]/255;
            otherwise
                error('unknown color');
        end
    else
        if nargin<4
            cmName = 'plot.viridis';
        elseif nargin<3
            error('not enough input arguments');
        elseif nargin>4
            error('to many input arguments');
        end
        cm = colormap(cmName);
        Pcm = linspace(0,1,numel(cm(:,1)))';
        if maxValue<minValue
            temp = minValue;
            minValue = maxValue;
            maxValue = temp;
        end
        if minValue==maxValue
            if value>minValue
                P = 1;
            elseif value<minValue
                P = 0;
            else
                P = 0.5;
            end
        else
            P = (value-minValue)/(maxValue-minValue);
            P = max([zeros(size(P)),min([ones(size(P)),P],[],2)],[],2);
        end
        rgb = interp1(Pcm,cm,P);
    end
end