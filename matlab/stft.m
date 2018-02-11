function [Y, F] = STFT(X, NFFT, WINDOW, NOVERLAP, Fs);

%
% Short-time Fourier Transform Analysis
%
%   [Y, F] = STFT(X, NFFT, WINDOW, NOVERLAP, Fs)
%
%     Y : stFT doamin output matrix (NFFT x L)
%     F : real frequencies corresponding to each FFT points
%     Fs : sampling freguency
%     NFFT : # of FFT points
%     WINDOW : window function
%     NOVERLAP : # of overlaped samples
%
%                      by Taesu Kim
%                           2003. 1. 26.
%

if nargin<5, Fs = 8000; end

WLEN = size(WINDOW, 1);
SHIFT = WLEN - NOVERLAP;

N = size(X, 1);
X = [zeros(WLEN, 1); X; zeros(NFFT, 1)];
%L = fix(N/SHIFT)+fix(WLEN/SHIFT);
L = fix((N+WLEN)/SHIFT) - 1;

Y = zeros(NFFT, L);

Xn = zeros(NFFT, 1);

for i = 1:L,
	sp = SHIFT*i + 1;
	Xn(1:WLEN) = WINDOW.*X(sp:sp+WLEN-1);
	Y(:, i) = fft(Xn);
end
Y = Y(1:fix(NFFT/2)+1,:);
% Y = Y(2:NFFT/2,:);
F = [0:fix(NFFT/2)]'.*Fs/NFFT;
