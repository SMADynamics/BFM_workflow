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
               rm_drift=False, 
               rm_drift_mode='linear', 
               rm_drift_pts=100, 
               rm_drift_plots=False, 
               stretch_xy=False, 
               stretch_xy_plots=False, 
               dist_beadsurf_wall=np.inf,
               filter_name='savgol', 
               filter_win=101,
               plots=False):
        ''' From x,y to torque.
            Take the (elliptical and drifting) rotating bead trajetory x(t),y(t), 
                1) correct drift
                2) rotate and stretch the trajectory into a circle
                3) calculate the angle Vs time
                4) calculate the angular speed Vs time 
                5) calculate the bead drag coefficient, with corrections due to surface proximity
                6) calculate and filter the torque trace
                7) plot everything for visual inspection
            
            Parameters:
            
            x, y                    : xy positions of the bead [pixel]
            bead_diam_m [1e-6]      : bead diameter
            FPS [1]                 : sampling rate, frames (points) per second. Used to calculate the speed.
            umppx [1]               : microns per pixel. Used to convert x,y to [m]eters. Use the defualt value of 1 if x,y are already given in [m].
            rm_drift [False]        : bool, remove drift by interpolation
            rm_drift_pts [100]      : number of pts to use in the interpolation for drift removal
            rm_drift_plots [False]  : plot the rm_drift traces
            rm_drift_mode           : 'linear' | 'spline'
            stretch_xy [False]      : bool, stretch the trajectory to a circle
            stretch_xy_plots        : plots relative to the xy stretch 
            dist_beadsurf_wall [inf]: gap between the surface of the bead and the wall, used to correct the drag
            filter_name ['savgol']  : 'savgol' or 'median', low pass filter for output torque trace
            filter_win [101]        : number of pts in filter window (odd number)
            plots [False]           : plot a summary of the results 
    
            Example:
            torque_pNnm_filtered = xy2torque.main(x, y, umppx=0.1475, bead_diam_m=1e-6, FPS=10000, filter_win=801, dist_beadsurf_wall=10e-9, rm_drift=True, rm_drift_plots=1, stretch_xy=1, stretch_xy_plots=1, plots=1)
    
        '''
        print('aaaaaaaaaaaaaaaaaaaaaaaaaaa')
        self.x=x
        self.y=y
        self.bead_diam_m=bead_diam_m,
        self.FPS=FPS, 
        self.umppx=umppx,
        self.rm_drift=rm_drift, 
        self.rm_drift_mode=rm_drift_mode, 
        self.rm_drift_pts=rm_drift_pts, 
        self.rm_drift_plots=rm_drift_plots
        self.stretch_xy=stretch_xy,  
        self.dist_beadsurf_wall=dist_beadsurf_wall,
        self.filter_name=filter_name, 
        self.filter_win=filter_win
        assert len(x)==len(y), 'Error: len(x) not equal to len(y)'
        print(self.rm_drift_pts)
        self.workflow()


    def workflow(self):
        # x,y to meters:
        x,y = self.x*self.umppx*1e-6, self.y*self.umppx*1e-6
    
        # 1) remove drift from x(t) and y(t):
        if self.rm_drift :
            print(f'main(): removing drift {self.rm_drift_mode}')
            x_rmdrift = filters.rm_interpolate(x, pts=self.rm_drift_pts, mode=self.rm_drift_mode, plots=self.rm_drift_plots, plot_signame='x')
            y_rmdrift = filters.rm_interpolate(y, pts=self.rm_drift_pts, mode=self.rm_drift_mode, plots=self.rm_drift_plots, plot_signame='y')
        else:
            print('main(): drift NOT removed.')
            x_rmdrift = x
            y_rmdrift = y
        
        # 2) stretch the elliptical trajectory (x(t),y(t)) into a circle:
        if self.stretch_xy:
            print('main(): stretching to circle')
            x_circ, y_circ = ellipse_fit.stretch_ellipse(x_rmdrift, y_rmdrift, plots=stretch_xy_plots)
        else:
            print('main(): NOT stretching to circle')
            x_circ = x_rmdrift
            y_circ = y_rmdrift
        
        # 3) calculate angle(t) [turns]:
        self.angle_turns = np.unwrap(np.arctan2(y - np.mean(y), x - np.mean(x)))/(2*np.pi)
        
        # 4) calculate speed [Hz]:
        self.speed_Hz = np.diff(self.angle_turns)*self.FPS
            
        # calc radius of trajectory [m]:
        self.traj_radius_m = np.hypot(x_circ, y_circ)  
        
        # 5) calculate the bead drag coefficient:
        self.drag_pNnms = drag.cal_drag(self.speed_Hz, self.traj_radius_m, self.bead_diam_m, self.dist_beadsurf_wall)
        print(f'main(): drag={self.drag_pNnms}')
    
        # 6) calculate the torque trace:
        self.torque_pNnm = self.drag_pNnms[:-1] * self.speed_Hz*2*np.pi
        print('main(): filtering torque trace')
        if self.filter_name == 'savgol':
            self.torque_pNnm_filter = filters.savgol_filter(self.torque_pNnm, self.filter_win, 5, mode=None)
        elif self.filter_name == 'median':
            self.torque_pNnm_filter = filters.median_filter(self.torque_pNnm, self.filter_win, fs=self.FPS )
        else: 
            self.torque_pNnm_filter = self.torque_pNnm
    
        # 7) plot everything:
        if plots:
            fig = plt.figure('xy2torque', clear=True)
            ax1 = fig.add_subplot(231)
            ax2 = fig.add_subplot(232)
            ax3 = fig.add_subplot(233)
            ax4 = fig.add_subplot(234)
            ax5 = fig.add_subplot(235)
            ax6 = fig.add_subplot(236)
            dw = 50 if len(x) > 1_000_000 else 1
            t = np.arange(len(x[::dw]))/FPS
            ax1.plot(x[::dw]*1e9, y[::dw]*1e9, ',', alpha=0.5, label='original')
            ax1.legend(fontsize=9)
            ax1.set_xlabel('x (nm)')
            ax1.set_ylabel('y (nm)')
            ax1.axis('equal')
            ax1.grid(0)
            ax2.plot(x_rmdrift[::dw]*1e9, y_rmdrift[::dw]*1e9, ',', alpha=0.5, label='drift corr.')
            ax2.legend(fontsize=9)
            ax2.axis('equal')
            ax2.set_xlabel('x (nm)')
            ax2.set_ylabel('y (nm)')
            ax2.grid(0)
            ax3.plot(x_circ[::dw]*1e9, y_circ[::dw]*1e9, ',', alpha=0.5, label='circ.stretched')
            ax3.legend(fontsize=9)
            ax3.axis('equal')
            ax3.set_xlabel('x (nm)')
            ax3.set_ylabel('y (nm)')
            ax3.grid(0)
            ax4.plot(t, self.angle_turns[::dw])
            ax4.set_ylabel('Angle (turns)')
            ax4.set_xlabel('Time (s)')
            ax4.grid(0)
            ax5.plot(t, self.traj_radius_m[::dw], '.', ms=2, alpha=0.2)
            ax5.set_ylabel('Traj. radius (m)')
            ax5.set_xlabel('Time (s)')
            ax5.grid(0)
            ax6.plot(t[:-1], self.torque_pNnm_filter[::dw], label=f'{filter_name} {filter_win}pts')
            ax6.legend(fontsize=9)
            ax6.set_ylabel('Torque (pN nm)')
            ax6.set_xlabel('Time (s)')
            ax6.grid(0)
            fig.tight_layout()
    
        print('main(): Done.')
        return torque_pNnm_filter
    
    
    
    
    
    
