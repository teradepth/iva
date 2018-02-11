# IVA: Independent Vector Analysis

## matlab implementations
ivabss.m
 
 Natural Gradient algorithm for Frecuency Domain Blind source separation based on Independent Vector Analysis
        
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

## python implementation

TO-DO


## References
[1] Taesu Kim, "Independent Vector Analysis" Ph.D. Dissertation, KAIST, 2007

[2] Taesu Kim, Hagai Attias, Soo-Young Lee, Te-Won Lee, "Blind source separation exploiting higher-order frequency dependencies" IEEE Transactions on Audio, Speech, and Language Processing 15 (1), 2007

[3] Intae Lee, Taesu Kim, Te-Won Lee, "Fast fixed-point independent vector analysis algorithms for convolutive blind source separation" Signal Processing 87 (8), 2007

[4] Taesu Kim, Torbjørn Eltoft, Te-Won Lee, "Independent vector analysis: An extension of ICA to multivariate components" International Conference on Independent Component Analysis and Signal Separation, 2006

[5] Taesu Kim, Intae Lee, Te-Won Lee, "Independent vector analysis: definition and algorithms", Fortieth Asilomar Conference on Signals, Systems and Computers, 2006
