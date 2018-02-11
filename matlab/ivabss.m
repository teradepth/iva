function [y, W] = ivabss(x, nfft, maxiter, tol, eta, nsou)

%  Natural Gradient algorithm for Frecuency Domain Blind source separation
%        based on Independent Vector Analysis
%
%    [y, W] = ivabss(x, nfft, maxiter, tol, eta, nsou)
%     y : separated signals (nsou x N)
%     W : unmixing matrices (nsou x nmic x nfft/2+1)
%     x : observation signals (nmic x N),
%           where nsou is # of sources, nmic is # of mics, and N is # of time frames
%     nfft : # of fft points (default =1024)
%     eta : learning rate (default =0.1)
%     maxiter : # of iterations (default =1000)
%     tol : When the difference of objective is less than tol,
%               the algorithm terminates (default =1e-6)
%     nsou : # of sources (default =nmic)
%
%
%                                               by Taesu Kim
%                                     Recently revised at Jan. 25, 2005

[nmic, nn] = size(x);
if ~exist('nfft','var')|isempty(nfft), nfft = 1024; end
if ~exist('eta','var')|isempty(eta), eta = 0.1; end
if ~exist('maxiter','var')|isempty(maxiter), maxiter = 1000; end
if ~exist('tol','var')|isempty(tol), tol=1e-6; end
if ~exist('nsou','var')|isempty(nsou), nsou = nmic; end

win = 2*hanning(nfft,'periodic')/nfft;
nol = fix(3*nfft/4);

for l=1:nmic,
    X(l,:,:) = conj(stft(x(l,:)', nfft, win, nol)');
end
clear x;

N = size(X,2);
nfreq = size(X,3);
epsi = 1e-6;
pObj = Inf;

% Meomory allocations
Wp = zeros(nsou,nsou,nfreq);
dWp = zeros(size(Wp));
Q = zeros(nsou,nmic,nfreq);
Xp = zeros(nsou,N,nfreq);
S = zeros(nsou,N,nfreq);
Ssq = zeros(nsou,N);
Ssq1 = zeros(nsou,N);

% Execute PCA and initialize
for k=1:nfreq,
    Xmean = mean(X(:,:,k),2)*ones(1,N);
    Rxx = (X(:,:,k)-Xmean)*(X(:,:,k)-Xmean)'/N;
    [E, D] = eig(Rxx);
    d = real(diag(D));
    [tmp, order] = sort(-d);
    E = E(:,order(1:nsou));
    D = diag(real(d(order(1:nsou)).^-.5));
    Q(:,:,k) = D*E';
    
    Xp(:,:,k) = Q(:,:,k)*(X(:,:,k)-Xmean);
    
    Wp(:,:,k) = eye(nsou);
end

% Start iterative learning algorithm
for iter=1:maxiter,

    dlw = 0;
    for k=1:nfreq,
        S(:,:,k) = Wp(:,:,k)*Xp(:,:,k);
    end
    
    Ssq = sum(abs(S).^2,3).^.5;
    Ssq1 = (Ssq+epsi).^-1;

    for k=1:nfreq,
        % Calculate multivariate score function and gradients
        Phi = Ssq1.*S(:,:,k);
        dWp(:,:,k) = (eye(nsou) - Phi*S(:,:,k)'/N)*Wp(:,:,k);
        dlw = dlw + log(abs(det(Wp(:,:,k)))+epsi);

    end
    
    % Update unmixing matrices
    Wp = Wp + eta*dWp;

    Obj = (sum(sum(Ssq))/N-dlw)/(nsou*nfreq);
    dObj = pObj-Obj;
    pObj = Obj;
    
    if mod(iter,10) == 0,
        fprintf('%d iterations: Objective=%e, dObj=%e\n',iter,Obj,dObj);
    end
    
    if abs(dObj)/abs(Obj) < tol, break; end
    
end

% Correct scaling of unmixing filter coefficients
for k=1:nfreq,
    W(:,:,k) = Wp(:,:,k)*Q(:,:,k);
    W(:,:,k) = diag(diag(pinv(W(:,:,k))))*W(:,:,k);
end

% % Spectral smoothing
% W(:,:,1) = (2*W(:,:,1) + W(:,:,2))/4;
% for k=2:nfreq-1,
%     W(:,:,k) = (W(:,:,k-1) + 2*W(:,:,k) + W(:,:,k+1))/4;
% end
% W(:,:,nfreq) = (W(:,:,nfreq-1) + 2*W(:,:,nfreq))/4;

%Calculate outputs
for k=1:nfreq,
    S(:,:,k) = W(:,:,k)*X(:,:,k);
end

% Re-synthesize the obtained source signals
for k=1:nsou,
    y(k,:) = istft(conj(squeeze(S(k,:,:))'), nn, win, nol)';
end
