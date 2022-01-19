%% path continuation - homotopy.f2
%
%   Institute of Dynamics and Vibration Research
%   Leibniz University Hannover
%   04.06.2020 - Alwin Förster
%
function [varargout] = f2(R,x,l,Opt)
    if Opt.jacobian
        [Rx,Jx] = R(x);
    else
        Rx = R(x);
    end
    f = 1+Rx.*Rx-l;
    varargout{1} = f;
    if Opt.jacobian
        nx = numel(x);
        J = [2*(Rx*ones(1,nx)).*Jx(1:nx,1:nx),-ones(nx,1)];
        varargout{2} = J;
    end
end