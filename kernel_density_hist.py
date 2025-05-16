import numpy as np
import matplotlib.pyplot as plt
from sklearn.neighbors import KernelDensity

def kernel_density_histo_old(sig, kernel='gaussian', logscale=False, band=1, return_all=False, plots=False):                                                                         
        ''' kernel density histogram of input sig.
        kernel: ['gaussian'|'tophat'|'epanechnikov'|'exponential'|'linear'|'cosine']
        band  : kernel bandwidth (~ bin size) 
        return (x, density) if return_all = False
        return (x, density, score, kde) if return_all = True
        '''
        if logscale:
            # TODO ok?
            Xout = np.logspace(np.log10(np.min(sig)*0.9), np.log10(np.max(sig)*1.1), 1000)[:, np.newaxis]
        else:
            Xout = np.linspace(np.min(sig)*0.9, np.max(sig)*1.1, 1000)[:, np.newaxis]
        kde = KernelDensity(kernel=kernel, bandwidth=band).fit(sig[:, np.newaxis])
        dens = np.exp(kde.score_samples(Xout))
        score = kde.score_samples
        if plots:
            plt.figure('kernel_density_histo', clear=True)
            plt.plot(Xout, dens)
        if return_all:
            return Xout, dens, score, kde
        else: 
            return Xout, dens


import numpy as np
import matplotlib.pyplot as plt
from sklearn.neighbors import KernelDensity

def kernel_density_histo(sig, kernel='gaussian', logscale=False, band=1.0, npts=1000, return_all=False, plots=False):
    """
    Compute and optionally plot the kernel density estimate (KDE) of a 1D signal.

    Parameters
    ----------
    sig : array_like
        1D input data.
    kernel : str, optional
        Kernel to use in the KDE. One of:
        ['gaussian', 'tophat', 'epanechnikov', 'exponential', 'linear', 'cosine'].
    logscale : bool, optional
        If True, x-axis is log-scaled. Only for positive-valued data.
    band : float, optional
        Bandwidth of the KDE (controls smoothing).
    return_all : bool, optional
        If True, return additional KDE objects and score method.
    plots : bool, optional
        If True, plot the KDE.

    Returns
    -------
    Xout : ndarray
        Evaluation points (shape: [N, 1]).
    dens : ndarray
        KDE evaluated at Xout.
    score : callable, optional
        KDE scoring function (only if return_all=True).
    kde : KernelDensity, optional
        Fitted KDE object (only if return_all=True).
    """
    sig = np.asarray(sig).ravel()
    if sig.size == 0:
        raise ValueError("Input `sig` is empty.")
    if logscale:
        if np.any(sig <= 0):
            raise ValueError("All signal values must be positive for logscale.")
        Xout = np.logspace(np.log10(np.min(sig) * 0.9), np.log10(np.max(sig) * 1.1), npts)[:, None]
    else:
        Xout = np.linspace(np.min(sig) * 0.9, np.max(sig) * 1.1, npts)[:, None]
    kde = KernelDensity(kernel=kernel, bandwidth=band)
    kde.fit(sig[:, None])
    log_dens = kde.score_samples(Xout)
    dens = np.exp(log_dens)
    if plots:
        plt.figure('kernel_density_histo', clear=True)
        plt.plot(Xout, dens, label=f'KDE ({kernel}, bw={band})')
        plt.xlabel('Value (log scale)' if logscale else 'Value')
        plt.ylabel('Density')
        plt.title('Kernel Density Estimate')
        plt.grid(True)
        plt.legend()
    if return_all:
        return Xout, dens, kde.score_samples, kde
    else:
        return Xout, dens

