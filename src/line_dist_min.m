%% path continuation - line_dist_min
%
%   Institute of Dynamics and Vibration Research
%   Leibniz University Hannover
%   26.11.2020 - Tido Kubatschek
% 
%   SRC: http://geomalgorithms.com/a07-_distance.html#dist3D_Segment_to_Segment
%
function [dist_min, tc, sc] = line_dist_min(p_1, p_2, q_1, q_2)

    if nargin ~= 4
        error('Not enough input arguments!');
    end

    u = p_2 - p_1;
    v = q_2 - q_1;
    w0 = p_1 - q_1;

    a = u'*u;
    b = u'*v;
    c = v'*v;
    d = u'*w0;
    e = v'*w0;

    sc = (b*e - c*d)/(a*c - b^2);
    tc = (a*e - b*d)/(a*c - b^2);
    tc = min(1, tc);
    sc = min(1, sc);
    dist_min = norm(w0 + (sc*u - tc*v));

    p_c = p_1 + sc*u;
    q_c = q_1 + tc*v;


end