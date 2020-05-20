%% path continuation - homotpy_merge
%
%   Institute of Dynamics and Vibration Research
%   Leibniz University Hannover
%   19.05.2020 - Alwin Förster
%
function [varargout] = homotopy_merge(G,R,x,l,Opt)
    if Opt.jacobian
        [Gx,JGx] = G(x);
        [Rx,JRx] = R(x);
        varargout{1} = (1-l)*Gx+l*Rx;
        varargout{2} = [(1-l)*JGx+l*JRx,-Gx+Rx];
    else
        varargout{1} = (1-l)*G(x)+l*R(x);
    end
end