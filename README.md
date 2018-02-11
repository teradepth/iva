# IVA
Independent Vector Analysis implementation

ivabss.m
  Natural Gradient algorithm for Frecuency Domain Blind source separation
        based on Independent Vector Analysis
    [y, W] = ivabss(x, nfft, maxiter, tol, eta, nsou)
     y : separated signals (nsou x N)
     W : unmixing matrices (nsou x nmic x nfft/2+1)
     x : observation signals (nmic x N),
           where nsou is # of sources, nmic is # of mics, and N is # of time frames
     nfft : # of fft points (default =1024)
     eta : learning rate (default =0.1)
     maxiter : # of iterations (default =1000)
     tol : When the difference of objective is less than tol,
               the algorithm terminates (default =1e-6)
     nsou : # of sources (default =nmic)

fiva.m
  Fast algorithm for Frecuency Domain Blind source separation
        based on Independent Vector Analysis

    [y, W] = fivabss(x, nfft, maxiter, tol, nsou)
     y : separated signals (nsou x N)
     W : unmixing matrices (nsou x nmic x nfft/2+1)
     x : observation signals (nmic x N),
           where nsou is # of sources, nmic is # of mics, and N is # of time frames
     nfft : # of fft points (default =1024)
     maxiter : # of iterations (default =1000)
     tol : When the increment of likelihood is less than tol,
               the algorithm terminates (deault =1e-6)
     nsou : # of sources (default =nmic)

