import numpy as np
import matplotlib.pyplot as plt

import filters
import drag
import ellipse_fit


class XY_2_Torque():

    def __init__(self, 
                 x, 
                 y, 
                 bead_diam_m=1e-6,
                 FPS=1, 
                 umppx=1,
                 correction_functions_order=['rm_outliers', 'rm_drift', 'stretch_xy'],
                 #rm_outliers=False,
                 rm_outliers_findparam=3,
                 rm_outliers_win=5,
                 rm_outliers_plots=False,
                 #rm_drift=False, 
                 rm_drift_mode='linear', 
                 rm_drift_pts=100, 
                 rm_drift_plots=False, 
                 #stretch_xy=False, 
                 stretch_xy_plots=False, 
                 dist_beadsurf_wall=50e-6,
                 filter_name='savgol', 
                 filter_win=101,
                 plots=False):
        ''' From x,y to torque. Take the (elliptical and drifting) rotating bead trajetory x(t),y(t) and:
                o) correct outlier points
                o) correct drift
                o) rotate and stretch the trajectory into a circle
                o) calculate the angle Vs time
                o) calculate the angular speed Vs time 
                o) calculate the bead drag coefficient, with corrections due to surface proximity
                o) calculate and filter the torque trace
                o) plot everything for visual inspection
            
            - Parameters:
            x, y                    : xy positions of the bead [pixel]
            bead_diam_m [1e-6]      : bead diameter
            FPS [1]                 : sampling rate, frames (points) per second. Used to calculate the speed.
            umppx [1]               : microns per pixel. Used to convert x,y to [m]eters. Use the defualt value of 1 if x,y are already given in [m].
            correction_functions_order ['rm_outliers', 'rm_drift', 'stretch_xy'] : list of correction functions names (if any), giving their order of execution
            #rm_outliers [False]       : remove oulier points in x,y
            rm_outliers_findparam [3]  : parameter to find outliers
            rm_outliers_win [3]        : window in pts to correct outliers        
            #rm_drift [False]          : bool, remove drift by interpolation
            rm_drift_pts [100]         : number of pts to use in the interpolation for drift removal
            rm_drift_plots [False]     : plot the rm_drift traces
            rm_drift_mode              : 'linear' | 'spline'
            #stretch_xy [False]        : bool, stretch the trajectory to a circle
            stretch_xy_plots           : plots relative to the xy stretch 
            dist_beadsurf_wall [50e-6] : gap between the surface of the bead and the wall, used to correct the drag [m]
            filter_name ['savgol']     : 'savgol' or 'median', low pass filter for output torque trace
            filter_win [101]           : number of pts in filter window (odd number)
            plots [False]              : plot a summary of the results 
    
            TODO
            - Example:
            torque_pNnm_filtered = xy2torque.main(x, y, umppx=0.1475, bead_diam_m=1e-6, FPS=10000, filter_win=801, dist_beadsurf_wall=10e-9, rm_drift=True, rm_drift_plots=1, stretch_xy=1, stretch_xy_plots=1, plots=1) 
        '''
        self.x=x
        self.y=y
        self.bead_diam_m=bead_diam_m
        self.FPS=FPS 
        self.umppx=umppx
        self.correction_functions_order=correction_functions_order
        #self.rm_outliers=rm_outliers
        self.rm_outliers_findparam=rm_outliers_findparam
        self.rm_outliers_win=rm_outliers_win
        self.rm_outliers_plots=rm_outliers_plots
        #self.rm_drift=rm_drift 
        self.rm_drift_mode=rm_drift_mode 
        self.rm_drift_pts=rm_drift_pts 
        self.rm_drift_plots=rm_drift_plots
        #self.stretch_xy=stretch_xy  
        self.stretch_xy_plots=stretch_xy_plots
        self.dist_beadsurf_wall=dist_beadsurf_wall
        self.filter_name=filter_name
        self.filter_win=filter_win
        self.plots=plots
        
        assert len(x)==len(y), 'Error: len(x) not equal to len(y)'
        self.correction_functions = {'rm_drift': self.remove_drift_funct, 'rm_outliers': self.remove_outliers_funct, 'stretch_xy': self.stretch_xy_funct}
        self.workflow()



    def workflow(self):
        # x,y to [meters]:
        x,y = self.x*self.umppx*1e-6, self.y*self.umppx*1e-6
        self.x_orig, self.y_orig = x, y
        
        # make corrections:
        self.x_rmdrift, self.y_rmdrift = x.copy(), y.copy()
        self.x_rmoutliers, self.y_rmoutliers = x.copy(), y.copy()  
        self.x_circ, self.y_circ = x.copy(), y.copy()
        for funct in self.correction_functions_order:
            x,y = self.correction_functions[funct](x,y)
        
        # calculate angle(t) [turns]:
        self.angle_turns = np.unwrap(np.arctan2(self.y_circ - np.mean(self.y_circ), self.x_circ - np.mean(self.x_circ)))/(2*np.pi)
        
        # calculate speed [Hz]:
        self.speed_Hz = np.diff(self.angle_turns)*self.FPS
            
        # calc radius of trajectory [m]:
        self.traj_radius_m = np.hypot(self.x_circ, self.y_circ)  
        
        # calculate the bead drag coefficient:
        self.drag_pNnms = drag.cal_drag(self.speed_Hz, self.traj_radius_m, self.bead_diam_m, self.dist_beadsurf_wall)
        print(f'XY_2_Torque.workflow(): Done angle(t), speed(t) traj_radius(t), drag(t) ({np.mean(self.drag_pNnms)} pNnms)')
    
        # calculate the torque trace:
        self.torque_pNnm = self.drag_pNnms[:-1] * self.speed_Hz*2*np.pi
        print(f'XY_2_Torque.workflow(): filtering torque and speed traces by {self.filter_name} win:{self.filter_win}')
        if self.filter_name == 'savgol':
            self.torque_pNnm_filter = filters.savgol_filter(self.torque_pNnm, self.filter_win, 5, mode=None)
            self.speed_Hz_filter = filters.savgol_filter(self.speed_Hz, self.filter_win, 5, mode=None)
        elif self.filter_name == 'median':
            self.torque_pNnm_filter = filters.median_filter(self.torque_pNnm, self.filter_win, fs=self.FPS )
            self.speed_Hz_filter = filters.median_filter(self.speed_Hz, self.filter_win, fs=self.FPS )
        elif self.filter_name == 'CK':
            self.torque_pNnm_filter = filters.CKfilter(self.torque_pNnm, K=5, M=self.filter_win, p=3, plots=False)
            self.speed_Hz_filter = filters.CKfilter(self.speed_Hz, K=5, M=self.filter_win, p=3, plots=False)
        else: 
            self.torque_pNnm_filter = self.torque_pNnm
            self.speed_Hz_filter = self.speed_Hz
    
        # plot everything:
        if self.plots:
            fig, d = plt.subplot_mosaic('''aabbcc\ndddeee\nffffhh\nggggjj''', num='xy2torque', clear=True)
            ax1 = d['a']
            ax2 = d['b']
            ax3 = d['c']
            ax4 = d['d']
            ax5 = d['e']
            ax6 = d['f']
            ax6a = d['h']
            ax7 = d['g']
            ax7a = d['j']
            t = np.arange(len(x))/self.FPS
            dw = 50 if len(x) > 1_000_000 else 1
            t = t[::dw]
            ax1.plot(self.x_orig[::dw]*1e9, self.y_orig[::dw]*1e9, ',', alpha=0.5, label='original')
            ax1.legend(fontsize=9)
            ax1.set_xlabel('x (nm)')
            ax1.set_ylabel('y (nm)')
            ax1.axis('equal')
            ax1.grid(0)
            ax2.plot(self.x_rmdrift[::dw]*1e9, self.y_rmdrift[::dw]*1e9, ',', alpha=0.5, label='drift corr.')
            ax2.legend(fontsize=9)
            ax2.axis('equal')
            ax2.set_xlabel('x (nm)')
            ax2.set_ylabel('y (nm)')
            ax2.grid(0)
            ax3.plot(self.x_circ[::dw]*1e9, self.y_circ[::dw]*1e9, ',', alpha=0.5, label='circ.stretched')
            ax3.legend(fontsize=9)
            ax3.axis('equal')
            ax3.set_xlabel('x (nm)')
            ax3.set_ylabel('y (nm)')
            ax3.grid(0)
            ax4.plot(t, self.angle_turns[::dw])
            ax4.set_ylabel('Angle (turns)')
            ax4.set_xlabel('Time (s)')
            ax4.grid(0)
            ax5.plot(t, self.traj_radius_m[::dw], ',', alpha=0.4)
            ax5.set_ylabel('Traj. radius (m)')
            ax5.set_xlabel('Time (s)')
            ax5.grid(0)
            try:
                ax6.plot(t[:-1], self.torque_pNnm_filter[::dw], ',', label=f'{self.filter_name} {self.filter_win}pts')
            except ValueError:
                ax6.plot(t, self.torque_pNnm_filter[::dw], ',', label=f'{self.filter_name} {self.filter_win}pts')
            ax6.legend(fontsize=8)
            ax6.set_ylabel('Torque (pN nm)')
            ax6.set_xlabel('Time (s)')
            ax6.grid(0)
            ax6a.hist(self.torque_pNnm_filter, 50, orientation='horizontal')
            ax6a.set_xticks([])
            ax6a.set_yticks([])
            try:
                ax7.plot(t[:-1], self.speed_Hz_filter[::dw], ',', label=f'{self.filter_name} {self.filter_win}pts')
            except ValueError:
                ax7.plot(t, self.speed_Hz_filter[::dw], ',', label=f'{self.filter_name} {self.filter_win}pts')
            ax7.legend(fontsize=8)
            ax7.set_ylabel('Speed (Hz)')
            ax7.set_xlabel('Time (s)')
            ax7.grid(0)
            ax7a.hist(self.speed_Hz_filter, 50, orientation='horizontal')
            ax7a.set_xticks([])
            ax7a.set_yticks([])

            fig.tight_layout()
        print('XY_2_Torque.workflow(): Done.')
    
    

    def remove_drift_funct(self, x, y):    
        print(f'remove_drift_funct(): removing drift {self.rm_drift_mode}')
        self.x_rmdrift = filters.rm_interpolate(x, pts=self.rm_drift_pts, mode=self.rm_drift_mode, plots=self.rm_drift_plots, plot_signame='x')
        self.y_rmdrift = filters.rm_interpolate(y, pts=self.rm_drift_pts, mode=self.rm_drift_mode, plots=self.rm_drift_plots, plot_signame='y')
        return self.x_rmdrift.copy(), self.y_rmdrift.copy()
    

    def remove_outliers_funct(self, x, y):
        print(f'remove_outliers_funct(): removing outliers (win:{self.rm_outliers_win} findparam:{self.rm_outliers_findparam})')
        self.x_rmoutliers = filters.outlier_smoother(x, m=4, win=3, plots=self.rm_outliers_plots, figname='rm_outliers x')[0]    
        self.y_rmoutliers = filters.outlier_smoother(y, m=4, win=3, plots=self.rm_outliers_plots, figname='rm_outliers y')[0]  
        return self.x_rmoutliers.copy(), self.y_rmoutliers.copy()


    def stretch_xy_funct(self, x, y): 
        print('stretch_xy_funct(): stretching x,y to circle')
        self.x_circ, self.y_circ = ellipse_fit.stretch_ellipse(x, y, plots=self.stretch_xy_plots)
        return self.x_circ.copy(), self.y_circ.copy()


