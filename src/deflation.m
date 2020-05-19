function [ Rdef ] = deflation( R, Xs, X )
    %% Deflation
    %
    p = 2;
    sigma = 10^-15;
    I = eye(length(R(Xs)));
    G = (1/sqrt((X-Xs)'*(X-Xs))^p+sigma)*I;
    Rdef = G*R(X);
    %
end