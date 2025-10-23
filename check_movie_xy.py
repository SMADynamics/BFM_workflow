### FP 251013
### check in a measurement everything about cw ccw measurements, from tracker to xy

import numpy as np
import matplotlib.pyplot as plt

import openTDMS
from tracked2torque import tracked_2_torque as t2t



class CheckVideo():

    def __init__(self, video_file_tdms='/home/francesco/PHYSMECHBIO/DATA/Simon/BA_campy/250926/CL_250926_111534.tdms', 
                 c0_trace=None, c1_trace=None, 
                 c0_trace_part=None, c1_trace_part=None, 
                 c0_frames=0, c1_frames=10, c2_frames=1):
        self.video_file_tdms = video_file_tdms
        self.c0_trace = c0_trace
        self.c1_trace = c1_trace
        self.c0_trace_part = c0_trace_part
        self.c1_trace_part = c1_trace_part
        self.c0_frames = c0_frames
        self.c1_frames = c1_frames
        self.c2_frames = c2_frames
        self.load_tracked_data()
        self.load_video_FPS()


    def load_video_FPS(self):
        ''' load video and FPS in tdms file '''
        print('load_video_FPS(): loading video ...')
        self.video, self.FPS = openTDMS.openTdmsOneROI(self.video_file_tdms, prints=True)


    def load_tracked_data(self):
        ''' load tracked 2 torque data'''
        self.trckd_path = self.video_file_tdms.strip('.tdms')+'/'
        self.tt = t2t(trckd_file=self.trckd_path,
                      roi_num=0,
                      c0=self.c0_trace, 
                      c1=self.c1_trace,
                      plots_xy=0,
                      bead_diam_m=1e-6,
                      dist_beadsurf_wall_m=100e-6,
                      umppx=0.1,
                      correction_functions_order=['rm_outliers', 'rm_drift', 'stretch_xy'],
                      rm_outliers_findparam=3,
                      rm_outliers_win=5,
                      rm_outliers_plots=False,
                      rm_drift_mode='linear',
                      rm_drift_pts=100,
                      rm_drift_qty='median',
                      rm_drift_qty_funct=np.median,
                      rm_drift_plots=True,
                      stretch_xy_plots=True,
                      filter_name='median',
                      filter_win=101,
                      plots_torque=True,
                      store_corr=False )

        
    def show_frames(self, c0=None, c1=None, c2=None, crop=(None,None,None,None)):
        ''' show cropped frames [c0:c1:c2] in video_file_tdms '''
        if c0==None: c0 = self.c0_frames
        if c1==None: c1 = self.c1_frames
        if c2==None: c2 = self.c2_frames
        nframes = (c1 - c0)/c2
        nframes = len(range(c0,c1,c2))
        xrowscols = int(np.ceil(np.sqrt(nframes)))
        print(f'show_video_frames(): numb. of frames: {xrowscols**2}, [{c0}:{c1}:{c2}]')
        fig_frames, axs_frames = plt.subplots(xrowscols, xrowscols, squeeze=False, num='show_video_frames', clear=True)
        i = 0
        for axr_frames in axs_frames:
            for ax_frames in axr_frames:
                ax_frames.imshow(self.video[c0+i, crop[0]:crop[1], crop[2]:crop[3]])
                ax_frames.set_xticks([])
                ax_frames.set_yticks([])
                ax_frames.set_title(c0+i, fontsize=5)
                i += c2
        fig_frames.tight_layout()
    
    
    def xy_plots(self, a=None, b=None, crop=(None,None,None,None), thr=0):
        ''' Plot cropped frames, xy, speed angle highlighting the region [a:b]. If b is small is easy to see CW Vs CCW.
                a, b : part of the traces to highlight
                crop : crop frames
                thr : threshold to plot omega > thr different color from omega < -thr
        '''
        tt = self.tt
        if a==None: 
            a = self.c0_trace_part
        else:
            self.c0_trace_part = a
        if b==None: 
            b = self.c1_trace_part
        else:
            self.c1_trace_part = b
        self.show_frames(c0=a, c1=b, c2=self.c2_frames, crop=crop)
        cond = (-thr < tt.speed_Hz_filter) * (tt.speed_Hz_filter < thr)
        fig = plt.figure('xy_plots', clear=1)
        fig.suptitle(self.video_file_tdms.rsplit('/')[-1], fontsize=8)
        ax1 = fig.add_subplot(311)
        ax1.plot(tt.x_corr[:-1][tt.speed_Hz_filter > thr ], tt.y_corr[:-1][tt.speed_Hz_filter > thr ], ',', c='0.5', alpha=0.3, label=(r'dark: $\omega >0$'))
        ax1.plot(tt.x_corr[:-1][tt.speed_Hz_filter < -thr], tt.y_corr[:-1][tt.speed_Hz_filter < -thr], ',', c=(0.8,0,0), alpha=0.3, label=(r'red :$\omega <0$'))
        ax1.plot(tt.x_corr[:-1][cond], tt.y_corr[:-1][cond], ',', c='g', alpha=0.5, label=(r'green: $\omega~0$'))
        ax1.plot(tt.x_corr[a:b], tt.y_corr[a:b], '-ob', ms=4)
        ax1.plot(tt.x_corr[b:b+1], tt.y_corr[b:b+1], 'bs', ms=9, label='end pt')
        ax1.legend()
        ax1.axis('equal')

        ax2b = fig.add_subplot(345)
        ax2b.plot(tt.angle_turns, 'y,')
        ax2b.plot(np.nonzero(cond)[0], tt.angle_turns[:-1][cond], 'g.', ms=2)
        ax2b.plot(range(a,b), tt.angle_turns[a:b], 'bo-', ms=4)
        ax2b.plot(b-1, tt.angle_turns[b-1], 'bs', ms=9)
        ax2b.set_xlabel('idx')
        ax2b.set_ylabel('angle turns')
        
        ax2 = fig.add_subplot(346)
        ax2.plot(np.mod(tt.angle_turns,1), 'y,')
        ax2.plot(np.nonzero(cond)[0], np.mod(tt.angle_turns[:-1][cond], 1), 'g.', ms=2)
        ax2.set_xlabel('idx')
        ax2.set_ylabel('angle mod 1 turn')
        
        ax2a = fig.add_subplot(347)
        ax2a.hist(np.mod(tt.angle_turns[:-1][cond], 1), 100, color='g', orientation='horizontal')
        ax2a.set_xlabel('Pauses count')
        ax2a.set_ylabel('angle mod 1 turn')
        
        ax21 = fig.add_subplot(348)
        ax21.plot(range(a-10, b+10), tt.angle_turns[a-10:b+10], 'yo-', ms=2, lw=1)
        ax21.plot(range(a,b), tt.angle_turns[a:b], 'bo-', ms=4)
        ax21.plot(b-1, tt.angle_turns[b-1], 'bs', ms=9)
        ax21.set_xlabel('idx')
        ax21.set_ylabel('angle turns')

        ax3 = fig.add_subplot(325)
        ax3.plot(tt.speed_Hz_filter, 'y,')
        ax3.plot(np.nonzero(cond)[0], tt.speed_Hz_filter[cond], 'g.')
        ax3.plot(range(a,b), tt.speed_Hz_filter[a:b], 'bo-', ms=4)
        ax3.plot(b-1, tt.speed_Hz_filter[b-1], 'bs', ms=9)
        ax3.set_xlabel('idx')
        ax3.set_ylabel('speed filtered Hz')
        ax3 = fig.add_subplot(326)
        ax3.plot(range(a-10,b+10), tt.speed_Hz_filter[a-10:b+10], 'y-o', ms=2, lw=1)
        ax3.plot(range(a,b), tt.speed_Hz_filter[a:b], 'bo-', ms=4)
        ax3.plot(b-1, tt.speed_Hz_filter[b-1], 'bs', ms=9)
        ax3.set_xlabel('idx')
        ax3.set_ylabel('speed filtered Hz')
        fig.tight_layout()



