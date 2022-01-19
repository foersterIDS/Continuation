%% path continuation - homotopy.fixNt
%
%   Institute of Dynamics and Vibration Research
%   Leibniz University Hannover
%   19.01.2022 - Alwin Förster
%
function [varargout] = fixNt(R,x,Rx0,x0,Opt)
    if Opt.jacobian
        [Rx,Jx] = R(x);
        nx = numel(x);
        Jx = Jx(1:nx,1:nx);
        varargout{1} = (x-x0).*(Rx-Rx0);
        varargout{2} = diag(Rx-Rx0)+((x-x0)*ones(1,nx)).*Jx;
    else
        Rx = R(x);
        varargout{1} = (x-x0).*(Rx-Rx0);
    end
end