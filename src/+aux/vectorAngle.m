%% path continuation - aux.vectorAngle
%
%   Institute of Dynamics and Vibration Research
%   Leibniz University Hannover
%   08.10.2020 - Tido Kubatschek
%
function [angle] = vectorAngle(u,v)
    %
    %% calculate cos(angle) using the dot product
    % must be in between -1 and 1
    %
    CosAngle = max(min(dot(u,v)/(norm(u)*norm(v)),1),-1);
    %
    %% calculate angle using acos
    %
    angle = real(acos(CosAngle));
    %
end

