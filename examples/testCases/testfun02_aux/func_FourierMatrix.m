%% Function: Fouriermatrices for fast transformation between time and frequency domain 
% 05.05.2021
% Florian Jaeger, IDS
% -------------------------------------------------------------------------
% Input:
%           nh      =   number of harmonics
%           nfft    =   number of fft points
% -------------------------------------------------------------------------
% Output:
%           G     =   Fouriermatrix time to frequency 
%           H     =   Fouriermatrix frequency to time 
% -------------------------------------------------------------------------

%%
function [G, H] = func_FourierMatrix(nh,nfft)
    idx_nh = [1:2:2*nh 2:2:2*nh];                               % Index sort [c1 s1 c2 s2...]
    tau = 0:2*pi/nfft:2*pi-2*pi/nfft;                           % nondimensional periodic time
    tau_nh = kron((1:nh)',tau);                                 % tau for every hamonic numbers 
    FFT_tau = [cos(tau_nh);sin(tau_nh)];
    H([1 idx_nh+1],:) = [ones(1,nfft); FFT_tau];                % frequency to time x(t)=X'*H
    G([1 idx_nh+1],:) = 1/nfft*[ones(1,nfft);2*FFT_tau];        % time to frequency X=G*x(t)
end