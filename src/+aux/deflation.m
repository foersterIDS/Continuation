%% path continuation - aux.deflation
%
%   Institute of Dynamics and Vibration Research
%   Leibniz University Hannover
%   08.05.2020 - Alwin Förster
%
function [ varargout ] = deflation( R, x_deflation, x, Opt )
    %% Deflation
    %
    p = 2;
    sigma = 10^-15;
    G = (1/sqrt((x-x_deflation)'*(x-x_deflation))^p+sigma);
    if Opt.jacobian
        [Rx,Jx] = R(x);
        JG = -p*(x-x_deflation)'*((x-x_deflation)'*(x-x_deflation))^(-p/2-1);
        varargout{1} = Rx*G;
        varargout{2} = Rx*JG+Jx*G;
    else
        Rx = R(x);
        varargout{1} = G*Rx;
    end
    %
end