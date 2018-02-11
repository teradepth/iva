function [Y] = ISTFT(X, N, WINDOW, NOVERLAP);

% Short-time Fourier Transform Synthesis by Overlap-add method
%
%   Y = ISTFT(X, N, WINDOW, NOVERLAP)
%
%     Y : time domain output vector
%     X : input stFT matrix (NFFT x L)
%     N : the length of output samples
%     WINDOW : window function
%     NOVERLAP : # of overlaped samples
%
%                      by Taesu Kim
%                           2003. 1. 26.
%

NFFT = (size(X, 1)-1)*2;
% NFFT = (size(X, 1)+1)*2;

L = size(X, 2);

WLEN = size(WINDOW, 1);
SHIFT = WLEN - NOVERLAP;

W = zeros(N+WLEN+NFFT, 1);
X = [X;conj(X(end-1:-1:2,:))];
% X = [zeros(1,L);X;zeros(1,L);conj(X(end:-1:1,:))];

Y = zeros(N+WLEN+NFFT, 1);

for i = 1:L,
	sp = SHIFT*i + 1;
	tmp = real(ifft(X(:, i)));
	
	W(sp:sp+WLEN-1) = W(sp:sp+WLEN-1) + WINDOW.^2;
	Y(sp:sp+WLEN-1) = Y(sp:sp+WLEN-1) + WINDOW.*tmp(1:WLEN);
	
% 	W(sp:sp+WLEN-1) = W(sp:sp+WLEN-1) + WINDOW;
% 	Y(sp:sp+NFFT-1) = Y(sp:sp+NFFT-1) + real(ifft(X(:, i)));
% 	Y(sp:sp+WLEN-1) = Y(sp:sp+WLEN-1) + tmp(1:WLEN);
end
W(SHIFT*(L+1)+1:SHIFT*(L+1)+WLEN) = W(SHIFT*(L+1)+1:SHIFT*(L+1)+WLEN) + WINDOW.^2;
% W(SHIFT*(L+1)+1:SHIFT*(L+1)+WLEN) = W(SHIFT*(L+1)+1:SHIFT*(L+1)+WLEN) + WINDOW;

Y = Y(WLEN+1:WLEN+N)./W(WLEN+1:WLEN+N);

