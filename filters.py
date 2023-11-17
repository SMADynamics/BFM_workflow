import numpy as np
import matplotlib.pyplot as plt


def savgol_filter(x, win, polyorder, mode=None, plots=0):
    ''' Savitzky-Gola filter of x,
            mode = 'valid' crop in [win:-win], 
            mode = None does nothing
    from https://stackoverflow.com/questions/20618804/how-to-smooth-a-curve-in-the-right-way/26337730#26337730'''
    import scipy.signal as sig
    if polyorder >= win:
        polyorder = win-1
        print('savgol_filter(): bad polyorder fixed to win-1')
    if np.mod(win,2) == 0:
        win = win+1
        print('Warning savgol_filter, win must be odd: forced win = '+str(win))
    y = sig.savgol_filter(x, window_length=win, polyorder=polyorder)
    if mode == 'valid':
        y = y[win:-win]
        x = x[win:-win]
    if plots:
        plt.figure('savgol_filter()')
        plt.clf()
        plt.plot(x, '.')
        plt.plot(y)
    return y



def median_filter(x, win=10, fs=1, usemode='reflect', plots=False):
    ''' median filter from scipy.ndimage '''
    from scipy.ndimage import median_filter
    y = median_filter(x, size=win, mode=usemode)
    if plots:
        plt.figure('median_filter()')
        plt.plot(np.arange(len(x))/fs, x, 'o')
        plt.plot(np.arange(len(y))/fs, y, '-')
        plt.xlabel('Time (s)')
    return y



def rm_interpolate(sig, p0=None, p1=None, pts=10, mode='spline', plot_signame='', plots=False):
    ''' remove interpolation made of pts from sig[p0:p1]
        mode = 'spline' or 'linear'
        see also: rm_interpolate_xy()
    '''
    # crop sig:
    sig = sig[p0:p1]
    # pick up pts points along sig by avg:
    idxs = np.linspace(0, len(sig)-1, pts, endpoint=1).astype(int)
    dd = np.diff(idxs)[0]
    sigpts = np.array([np.median(sig[i:i + dd]) for i in idxs])
    sigpts[-1] = np.median(sig[-dd:])
    sigidx = range(0, len(sig))
    if mode == 'spline':
        from scipy.interpolate import splev
        from scipy.interpolate import splrep
        # interpolate by spline:
        spline = splrep(idxs, sigpts)
        interp = splev(sigidx, spline)
        # remove interpolation from sig :
        sig_out = sig - interp + sigpts[0]
    elif mode == 'linear':
        from scipy.signal import detrend
        sig_out = detrend(sig) + sigpts[0]
        interp = sig - sig_out + sigpts[0]
    if plots:
        if len(sig) > 1000000:
            dw = 50
        else:
            dw = 1
        fig = plt.figure('rm_interpolate '+plot_signame, clear=True)
        ax1 = fig.add_subplot(211)
        ax2 = fig.add_subplot(212)
        ax1.plot(sigidx[::dw], sig[::dw], '.', ms=2, label='raw data', alpha=0.2)
        ax1.plot(idxs, sigpts, 'o', mfc='none', label='interp. pts')
        ax1.plot(sigidx[::dw], interp[::dw], label='interpolation')
        ax1.set_ylabel('sig', fontsize=14)
        ax1.legend()
        ax2.plot(sigidx[::dw], sig_out[::dw], '.', ms=2, alpha=0.2, label='corrected')
        ax2.set_ylabel('sig', fontsize=14)
        ax2.set_xlabel('index', fontsize=14)
        ax2.legend()
    return sig_out



