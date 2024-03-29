%% path continuation - homotopy.merge
%
%   Institute of Dynamics and Vibration Research
%   Leibniz University Hannover
%   19.05.2020 - Alwin F�rster
%
function [varargout] = merge(G,R,x,l,Opt)
    if Opt.jacobian
        [Gx,JGx] = G(x);
        [Rx,JRx] = R(x);
        JGx = JGx(1:numel(Rx),1:numel(x));
        JRx = JRx(1:numel(Rx),1:numel(x));
        varargout{1} = (1-l)*Gx+l*Rx;
        varargout{2} = [(1-l)*JGx+l*JRx,-Gx+Rx];
    else
        varargout{1} = (1-l)*G(x)+l*R(x);
    end
end