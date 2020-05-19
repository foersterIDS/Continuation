function [ varargout ] = deflation( R, x_deflation, x )
    %% Deflation
    %
    p = 2;
    sigma = 10^-15;
    G = (1/sqrt((x-x_deflation)'*(x-x_deflation))^p+sigma);
    if abs(nargout(R))==2
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