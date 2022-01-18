%% path continuation - aux.vector_angle
%
%   Institute of Dynamics and Vibration Research
%   Leibniz University Hannover
%   08.10.2020 - Tido Kubatschek
%
function [angle] = vector_angle(u,v)
    %
    %% calculate cos(angle) using the dot product
    % must be in between -1 and 1
    %
    Cos_angle = max(min(dot(u,v)/(norm(u)*norm(v)),1),-1);
    %
    %% calculate angle using acos
    %
    angle = real(acos(Cos_angle));
    %
end

