%% path continuation - homotopy.newton
%
%   Institute of Dynamics and Vibration Research
%   Leibniz University Hannover
%   19.05.2020 - Alwin F�rster
%
function [varargout] = newton(R,x,x0)
    if abs(nargout(R))==2
        [Rx,Jx] = R(x);
        varargout{1} = Rx-R(x0);
        varargout{2} = Jx;
    else
        varargout{1} = R(x)-R(x0);
    end
end