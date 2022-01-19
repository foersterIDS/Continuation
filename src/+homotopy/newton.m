%% path continuation - homotopy.newton
%
%   Institute of Dynamics and Vibration Research
%   Leibniz University Hannover
%   19.05.2020 - Alwin Förster
%
function [varargout] = newton(R,x,Rx0,Opt)
    if Opt.jacobian
        [Rx,Jx] = R(x);
        Jx = Jx(1:numel(Rx),1:numel(x));
        varargout{1} = Rx-Rx0;
        varargout{2} = Jx;
    else
        varargout{1} = R(x)-Rx0;
    end
end