%% path continuation - aux.lineDistMin
%
%   Institute of Dynamics and Vibration Research
%   Leibniz University Hannover
%   26.11.2020 - Tido Kubatschek
% 
%   SRC: http://geomalgorithms.com/a07-_distance.html#dist3D_SegmentTo_Segment
%
function [distMin, tc, sc, pI] = lineDistMin(p1, p2, q1, q2)
    %
    if nargin ~= 4
        error('Not enough input arguments!');
    end
    %
    %% Define variables 
    %
    u = p2 - p1;
    v = q2 - q1;
    w0 = p1 - q1;
    %
    a = u'*u;
    b = u'*v;
    c = v'*v;
    d = u'*w0;
    e = v'*w0;
    %
    %% Calc min. distance
    %
    sc = (b*e - c*d)/(a*c - b^2);
    tc = (a*e - b*d)/(a*c - b^2);
    tc = min(1, tc);
    sc = min(1, sc);
    %
    distMin = norm(w0 + (sc*u - tc*v));
    %
    %% Calc point of intersection
    %
    pI = p1 + sc*u;
    %
end