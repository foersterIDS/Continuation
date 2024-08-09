%% path continuation - aux.deflation
%
%   Institute of Dynamics and Vibration Research
%   Leibniz University Hannover
%   08.05.2020 - Alwin Förster
%
function [ varargout ] = deflation( R, xDeflation, x, hasJacobian )
    %% Deflation
    %
    p = 2;
    sigma = 10^-15;
    G = (1/sqrt((x-xDeflation)'*(x-xDeflation))^p+sigma);
    if hasJacobian
        [Rx,Jx] = R(x);
        JG = -p*(x-xDeflation)'*((x-xDeflation)'*(x-xDeflation))^(-p/2-1);
        varargout{1} = Rx*G;
        varargout{2} = Rx*JG+Jx*G;
    else
        Rx = R(x);
        varargout{1} = G*Rx;
    end
    %
end