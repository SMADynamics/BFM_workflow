import numpy as np
import matplotlib.pyplot as plt


#TODO
def correct_sig_modulation(sig, angle_turns, method='poly', polydeg=10, add='mean', interp_pts=100, plots=False, plots_figname='', plot_ss=30, plots_test=False, return_all=False):
    '''correct the 1-turn periodic modulation (for a signal of the BFM eg: omega, z, radius), 
       by removing a polynome fit of sig Vs mod(angle_turns,1).
         sig, angle_turns: any signal and relative angle_turns trace (they have same length)
         method: 'poly' polynome fit, 'interp' interpolation
         polydeg : polyn degree to fit
         add : ['mean', 'max'] value to add to the corrected signal. If add==None, then np.mean(sig) is added. But if sig=radius of xy traj, it should  be different TODO what?
         interp_pts : numb of point to use to interpolate
         plot_ss: subsample when plotting
    return signal corrected
    '''
    if add == 'mean':
        _add = np.mean(sig)
    # angle turns in 0,1:
    am  = np.mod(angle_turns - angle_turns[0], 1)
    if method == 'poly':
        # polyn fit:
        pf = np.polyfit(am, sig, polydeg)
        po = np.poly1d(pf)
        if add == 'max':
            _add = np.max(po(am))
        # sig corrected:
        sig_corr = sig - po(am) + _add
    elif method == 'interp':
        sig_corr, interp_f = rm_interpolate_xy(am, sig, npts=interp_pts, add=add, plots=plots)
    if plots :
        if plots_figname:
            plt.figure(plots_figname, clear=True)
        else:
            plt.figure('correct_sig_modulation', clear=True)
        plt.subplot(321)
        plt.plot(sig[::plot_ss], '-', alpha=0.8, label='sig.raw')
        plt.plot(sig_corr[::plot_ss], '-', alpha=0.6, label='sig.corr')
        plt.legend(fontsize=9)
        plt.subplot(322)
        plt.plot(angle_turns[::plot_ss], '.', ms=2, label='angle_turns')
        plt.legend(fontsize=9)
        plt.subplot(312)
        if method == 'poly':
            plt.plot(am[::plot_ss], sig[::plot_ss], ',', ms=1, alpha=0.3, label='sig.raw')
            plt.plot(am[::plot_ss], po(am)[::plot_ss], ',', ms=2, label='poly.fit')
        plt.xlabel('angle_turns mod 1')
        plt.legend(fontsize=9)
        plt.subplot(313)
        if method == 'poly':
            plt.plot(am[::plot_ss], sig_corr[::plot_ss], ',', ms=1, alpha=0.3, label='sig.corr.')
        plt.xlabel('angle_turns mod 1')
        plt.legend(fontsize=9)
        plt.tight_layout()
    return sig_corr



def CKfilter(signal, K=5, M=5, p=3, plots=0):
    ''' by AN 2025
    This code was adapted from my Matlab code, 2012
    Algorithm from Chung & Kennedy,
    Forward-backward non-linear filtering technique for extracting small biological signals from noise,
    J Neurosci Meth, 1991
    M is size of the analysis window, p is a weighting factor that typically ranges
    between 1 and 100, and K is the number of forward and backward predictors
    '''
    box = np.concatenate((np.ones(int(M))/int(M), np.zeros(int(M-1))))
    bkwd_mn = np.convolve(signal, box, mode='same')
    for_mn = np.convolve(signal, box[::-1], mode='same')
    f_i, b_i = np.zeros((K,len(signal))), np.zeros((K,len(signal)))
    for i in range(1,K+1):
        new_box = box*M
        for_spread = (signal - for_mn)**2
        bkwd_spread = (signal - bkwd_mn)**2
        f_sum = np.convolve(for_spread,new_box[::-1],mode='same')
        b_sum = np.convolve(bkwd_spread,new_box,mode='same')
        f_i[i-1,:] =  1/(2*K)*(f_sum)**-p
        b_i[i-1,:] =  1/(2*K)*(b_sum)**-p
    #normalize f_i and b_i
    f_i, b_i = f_i/np.sum(f_i+b_i,axis=0), b_i/np.sum(f_i+b_i,axis=0)
    filtered_sig = np.zeros(len(for_mn))
    for i in range(K):
        filtered_sig = filtered_sig + f_i[i]*for_mn + b_i[i]*bkwd_mn
    if plots:
        plt.figure('CKfilter', clear=True)
        plt.plot(signal, '.', ms=1)
        plt.plot(filtered_sig)
    return filtered_sig



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
        if len(sig) > 500000:
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
        ax1.legend(fontsize=9)
        ax2.plot(sigidx[::dw], sig_out[::dw], '.', ms=2, alpha=0.2, label='corrected')
        ax2.set_ylabel('sig', fontsize=14)
        ax2.set_xlabel('index', fontsize=14)
        ax2.legend(fontsize=9)
    return sig_out



def outlier_smoother(x, m=3, win=3, plots=False, figname=''):
    ''' finds outliers in x (points > m*mdev(x)) [mdev:median deviation] 
    and replaces them with the median of win points around them 
    return x_corrected and number of outliers found '''
    x_corr = np.copy(x)
    d = np.abs(x - np.median(x))
    mdev = np.median(d)
    idxs_outliers = np.nonzero(d > m*mdev)[0]
    if len(idxs_outliers): print(f'outlier_smoother(): removing {len(idxs_outliers)} outliers [win:{win}]...')
    k = 0
    for i in idxs_outliers:
        if 100%100 == 0: print(f'outlier_smoother(): {k}/{len(idxs_outliers)}', end='\r')
        if i-win < 0:
            x_corr[i] = np.median(np.append(x[0:i], x[i+1:i+win+1]))
        elif i+win+1 > len(x):
            x_corr[i] = np.median(np.append(x[i-win:i], x[i+1:len(x)]))
        else:
            x_corr[i] = np.median(np.append(x[i-win:i], x[i+1:i+win+1]))
        k += 1
    if k>0: print()
    if plots:
        dw = 50 if len(x) > 500_000 else 1
        fig = plt.figure(f'outlier_smoother '+figname, clear=True)
        ax1 = fig.add_subplot(211)
        ax1.plot(np.arange(0, len(x), dw), x[::dw], label='orig.', lw=2)
        ax1.plot(idxs_outliers, x[idxs_outliers], 'ro', label='outliers')
        ax1.legend(loc='upper left', framealpha=0.3, fontsize=9)
        ax2 = fig.add_subplot(212, sharex=ax1)
        ax2.plot(np.arange(0, len(x), dw), x_corr[::dw], '-o', label='corrected')
        ax2.legend(loc='upper left', framealpha=0.3, fontsize=9)
    return x_corr, len(idxs_outliers)

