import numpy as np
import matplotlib.pyplot as plt
import scipy.signal


##################################################################################################################
## used for verapamil, and possibly new here in filters.py. 
# Remove pauses in xy by clickin in a gui or by differenet other ways:

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



def rm_regions_gui(x, plot_signame=''):
    ''' GUI to select regions to remove from x. Click on the plot to select the beginning and end of each region to remove. Return a list of tuples with the indexes of the regions to remove. '''
    fig = plt.figure('rm_regions_gui ' + plot_signame, clear=True)
    plt.plot(x)
    plt.title('2 clicks beg/end of each region to remove. Left: add. Right: rm. Middle: done.', fontsize=9)
    plt.xlabel('idx', fontsize=12)
    plt.ylabel('x', fontsize=12)
    plt.pause(0.1)
    print('rm_regions_gui(): Click on the plot to select regions to remove...')
    idxs = plt.ginput(n=-1, timeout=0)
    if len(idxs) == 0:
        print('rm_regions_gui(): no points selected, no regions to remove.')
        return []
    idxs = np.array(idxs).astype(int)[:,0]
    print(f'rm_regions_gui(): idxs: {idxs}')
    idxs_filtered = [(idxs[i], idxs[i+1]) for i in range(0, len(idxs)-1, 2)]
    # plot regions selected:
    if idxs_filtered:
        plt.figure('rm_regions_gui ' + plot_signame)
        for a,b in idxs_filtered:
            plt.axvspan(a, b, color='r', alpha=0.2, ls='-', lw=2)
    else:
        idxs_filtered = []
    plt.pause(0.1)
    return idxs_filtered



def rm_interpolate(sig, p0=None, p1=None, rm_regions=[], rm_regions_thr=0.5, rm_regions_extend=0, wins=10, mode='polyfit', polyfit_deg=3, qty='funct', qty_funct=np.median, plot_signame='', plots=False, plot_suptitle=None, return_metadata=False):
    ''' Remove from a signal its interpolation. 
        Divide sig in time windows ('wins'). In each, calculate a quantity (qty), eg: median, to produce points for the interpolation. Then interpolate these points and remove the interpolation from sig.
            p0,p1 : idx to crop sig[p0:p1] before interpolation.
            rm_regions : remove regions of sig between indexes a_i and b_i, before interpolation. Useful to remove pauses.
                         Either a list of tuples [(a1,b1),(a2,b2),...] to remove from sig before calculating the interpolation. 
                         Or str 'auto', to use find_pauses(..., thr=rm_regions_thr) to automatically find pauses.
                         Or str 'gui', to use rm_regions_gui() to manually select regions to remove.
                         rm_regions_extend : extend the regions to remove by this factor (eg: 0.5 means extend by 50% before and after the detected region).
            qty : ['median'|'min'|'max'|'funct'] the quantity to calculate in each window
                  to produce the points for the interpolation.
            qty_funct : if qty='funct', then 'qty_funct' must be a function (eg: qty_funct = np.mean)
            mode : ['spline' | 'linear' | 'polyfit' | 'none' ] method to use for interpolation. 'none' means no correction, just return the cropped sig.
            polyfit_deg : if mode=='polyfit', degree of the polynome to fit
            
        see also: 
        rm_interpolate_xy()

        TODO: the way an interpolating point is placed in its window can be improved, now it is placed at the beginning of the window, but it could be placed in the middle.
    '''
    # crop sig:
    sig = sig[p0:p1]
    # pick up pts points along sig:
    idxs = np.linspace(0, len(sig)-1, wins, endpoint=1).astype(int)
    dd = np.diff(idxs)[0]   
    # remove points in idxs that are in rm_regions:
    if isinstance(rm_regions, list):
        rm_regions_idxs = rm_regions
    elif rm_regions == 'auto':
        rm_regions_idxs = find_pauses(sig, thr=rm_regions_thr, extend=rm_regions_extend, plots=plots, NFFT=dd, plot_signame=plot_signame)
    elif rm_regions == 'gui':
        rm_regions_idxs = rm_regions_gui(sig, plot_signame=plot_signame)
    else:
        raise ValueError('rm_interpolate(): ERROR "rm_regions" not well defined.')
    for a, b in rm_regions_idxs:
        idxs = idxs[(idxs < a) | (idxs >= b)]
    print(f'rm_interpolate(): removed regions: {rm_regions_idxs}')
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
        raise ValueError('rm_interpolate(): ERROR "qty" not well defined.')
    sigidx = range(0, len(sig))
    if mode == 'spline':
        from scipy.interpolate import splev
        from scipy.interpolate import splrep
        # interpolate by spline:
        spline = splrep(idxs, sigpts, k=1)
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
    elif mode == 'none':
        sig_out = sig
        interp = np.zeros_like(sig)
    else:
        raise ValueError('rm_interpolate_xy(): ERROR "mode" not well defined.')
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
        plt.pause(0.1)
    if return_metadata:
        metadata = {'qty': qty, 'wins': wins, 'mode': mode, 'polyfit_deg': polyfit_deg, 'rm_regions_idxs': [(int(a), int(b)) for a, b in rm_regions_idxs]}
        return sig_out, metadata
    else:
        return sig_out



#######################################################################################################################
# speed modulation correction, possibly in filters.py :

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

