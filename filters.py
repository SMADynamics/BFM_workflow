import numpy as np
import matplotlib.pyplot as plt
import scipy.signal








#TODO adjust and use:
# def correct_sig_modulation(sig, angle_turns, method='poly', polydeg=10, add='mean', interp_pts=100, plots=False, plots_figname='', plot_ss=30, plots_test=False, return_all=False):
#     '''correct the 1-turn periodic modulation (for a signal of the BFM eg: omega, z, radius), 
#        by removing a polynome fit of sig Vs mod(angle_turns,1).
#          sig, angle_turns: any signal and relative angle_turns trace (they have same length)
#          method: 'poly' polynome fit, 'interp' interpolation
#          polydeg : polyn degree to fit
#          add : ['mean', 'max'] value to add to the corrected signal. If add==None, then np.mean(sig) is added. But if sig=radius of xy traj, it should  be different TODO what?
#          interp_pts : numb of point to use to interpolate
#          plot_ss: subsample when plotting
#     return signal corrected
#     '''
#     if add == 'mean':
#         _add = np.mean(sig)
#     # angle turns in 0,1:
#     am  = np.mod(angle_turns - angle_turns[0], 1)
#     if method == 'poly':
#         # polyn fit:
#         pf = np.polyfit(am, sig, polydeg)
#         po = np.poly1d(pf)
#         if add == 'max':
#             _add = np.max(po(am))
#         # sig corrected:
#         sig_corr = sig - po(am) + _add
#     elif method == 'interp':
#         sig_corr, interp_f = rm_interpolate_xy(am, sig, npts=interp_pts, add=add, plots=plots)
#     if plots :
#         if plots_figname:
#             plt.figure(plots_figname, clear=True)
#         else:
#             plt.figure('correct_sig_modulation', clear=True)
#         plt.subplot(321)
#         plt.plot(sig[::plot_ss], '-', alpha=0.8, label='sig.raw')
#         plt.plot(sig_corr[::plot_ss], '-', alpha=0.6, label='sig.corr')
#         plt.legend(fontsize=9)
#         plt.subplot(322)
#         plt.plot(angle_turns[::plot_ss], '.', ms=2, label='angle_turns')
#         plt.legend(fontsize=9)
#         plt.subplot(312)
#         if method == 'poly':
#             plt.plot(am[::plot_ss], sig[::plot_ss], ',', ms=1, alpha=0.3, label='sig.raw')
#             plt.plot(am[::plot_ss], po(am)[::plot_ss], ',', ms=2, label='poly.fit')
#         plt.xlabel('angle_turns mod 1')
#         plt.legend(fontsize=9)
#         plt.subplot(313)
#         if method == 'poly':
#             plt.plot(am[::plot_ss], sig_corr[::plot_ss], ',', ms=1, alpha=0.3, label='sig.corr.')
#         plt.xlabel('angle_turns mod 1')
#         plt.legend(fontsize=9)
#         plt.tight_layout()
#     return sig_corr


# TODO: careful I see arteftacts, at the beg and end of the filtered trace, few points out of the trace, going down. rm avg does not rm arteftacts.
def CKfilter(signal, K=5, M=5, p=3, plots=0):
    ''' by AN 2025. This code was adapted from my Matlab code, 2012
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



def find_pauses(x, thr=0.7, extend=0, NFFT=2**10, plots=False, plot_signame=''):
    ''' find pauses in sinusoidal-like signal x, by  thresholding the sum(axis=0) of the spectrogram'''
    f, t, Sxx = scipy.signal.spectrogram(x, nperseg=NFFT, noverlap=0)
    m = np.mean(Sxx, axis=0)
    m = m / np.max(m)
    idxs = np.nonzero(m < thr)[0]
    # find indexes in idxs that are not consecutive, and keep only the first and last of each group:
    idxs_groups = np.split(idxs, np.where(np.diff(idxs) > 1)[0] + 1)
    idxs_filtered = [(int(group[0]*(1 - extend)), int(group[-1]*(1 + extend))) for group in idxs_groups if len(group) > 0]
    sample_start = np.arange(m.shape[0]) * NFFT
    #print(f'find_pauses(): found idx_filtered: {type(idxs_filtered)} {idxs_filtered}')
    if idxs_filtered is not None:
        idxs_out = sample_start[idxs_filtered]
    else:
        idxs_out = np.array([])
    print(f'find_pauses(): found {len(idxs_out)} pauses')
    if plots:
        plt.figure('find_pauses()' + plot_signame, clear=True)
        plt.subplot(211)
        plt.plot(x)
        plt.plot(sample_start[idxs], x[sample_start[idxs]], 'ro', mfc='none', label='pauses')
        plt.ylabel('input signal', fontsize=12)
        plt.subplot(212)
        plt.plot(sample_start, m, label='mean spectrogram')
        plt.plot(sample_start[idxs], m[idxs], 'ro', mfc='none', label='pauses')
        plt.axhline(thr, color='k', ls='--', label='threshold')
        for a,b in idxs_filtered:
            plt.axvspan(sample_start[a], sample_start[b]+NFFT, color='r', alpha=0.2, ls='-', lw=2)
        plt.xlabel('idx', fontsize=12)
        plt.ylabel('mn spectrogram', fontsize=12)
        plt.legend(fontsize=9)
    return idxs_out



# WARNING added 'qty' and 'qty_funct'. TODO update callers:
def rm_interpolate(sig, p0=None, p1=None, rm_regions=[], rm_regions_thr=0.5, rm_regions_extend=0, wins=10, mode='polyfit', polyfit_deg=3, qty='funct', qty_funct=np.median, plot_signame='', plots=False, plot_suptitle=None):
    ''' Remove from a signal its interpolation. 
        Divide sig in time windows ('wins'). In each, calculate a quantity (qty), eg: median, to produce points for the interpolation. Then interpolate these points and remove the interpolation from sig.
            p0,p1 : idx to crop sig[p0:p1] before interpolation.
            rm_regions : remove regions of sig between indexes a_i and b_i, before interpolation. Useful to remove pauses.
                         Either a list of tuples [(a1,b1),(a2,b2),...] to remove from sig before calculating the interpolation. 
                         Or str 'auto', to use find_pauses(..., thr=rm_regions_thr) to automatically find pauses.
            qty : ['median'|'min'|'max'|'funct'] the quantity to calculate in each window
                  to produce the points for the interpolation.
            qty_funct : if qty='funct', then 'qty_funct' must be a function (eg: qty_funct = np.mean)
            mode : ['spline' | 'linear' | 'polyfit'] method to use for interpolation
            polyfit_deg : if mode=='polyfit', degree of the polynome to fit
            
        see also: 
        rm_interpolate_xy()

        TODO: the way a point in its window is placed can be improved, now it is placed at the beginning of the window, but it could be placed in the middle or at the end.
    '''
    # crop sig:
    sig = sig[p0:p1]
    # pick up pts points along sig:
    idxs = np.linspace(0, len(sig)-1, wins, endpoint=1).astype(int)
    dd = np.diff(idxs)[0]   
    # remove points in idxs that are in rm_regions:
    if rm_regions and isinstance(rm_regions, list):
        for a, b in rm_regions:
            idxs = idxs[(idxs < a) | (idxs >= b)]
    elif rm_regions == 'auto':
        idxs_pauses = find_pauses(sig, thr=rm_regions_thr, extend=rm_regions_extend, plots=plots, NFFT=dd, plot_signame=plot_signame)
        for a, b in idxs_pauses:
            idxs = idxs[(idxs < a) | (idxs >= b)]
    # calc. qty and interp. points in each window:
    if qty == 'median':
        sigpts = np.array([np.median(sig[i:i + dd]) for i in idxs])
        sigpts[-1] = np.median(sig[-dd:])
    elif qty == 'min':
        sigpts = np.array([np.min(sig[i:i + dd]) for i in idxs])
        sigpts[-1] = np.min(sig[-dd:])
    elif qty == 'max':
        sigpts = np.array([np.max(sig[i:i + dd]) for i in idxs])
        sigpts[-1] = np.max(sig[-dd:])
    elif qty == 'funct':
        sigpts = np.array([qty_funct(sig[i:i + dd]) for i in idxs])
        sigpts[-1] = qty_funct(sig[-dd:])
    else:
        print('rm_interpolate_xy(): ERROR "qty" not well defined.')
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
    elif mode == 'polyfit':
        pf = np.polyfit(idxs, sigpts, deg=polyfit_deg)
        po = np.poly1d(pf)
        interp = po(sigidx)
        sig_out = sig - interp + sigpts[0]
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
        ax1.legend(fontsize=9, labelspacing=0)
        ax2.plot(sigidx[::dw], sig_out[::dw], '.', ms=2, alpha=0.2, label='corrected')
        ax2.set_ylabel('sig', fontsize=14)
        ax2.set_xlabel('index', fontsize=14)
        ax2.legend(fontsize=9)
        fig.suptitle(plot_suptitle, fontsize=8)
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


if __name__ == '__main__':
    import tracked2torque
    tt = tracked2torque.tracked_2_torque(trckd_file='/home/francesco/CL_241202_144019/12/trajectory.pt', plots_torque=0)
    rm_interpolate(tt.x_orig, pts=5, rm_ab=[(80000,520000)], mode='spline', qty='median', plots=1)