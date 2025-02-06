# calls to get torque (and speed, filtered) from the tracked .tdms file (both FP or TF tracker)
import tracked2xy
import xy2torque



def tracked_2_torque(trckd_file='', 
                     roi_num=0,
                     plots_xy=False,
                     bead_diam_m=1e-6, 
                     dist_beadsurf_wall=100e-6,
                     umppx=0.1,
                     rm_drift=False, 
                     rm_drift_mode='linear',
                     rm_drift_pts=100, 
                     rm_drift_plots=False,
                     stretch_xy=False, 
                     stretch_xy_plots=False,
                     filter_name='median', 
                     filter_win=101,
                     plots_torque=False ):
    ''' calls to get torque (and speed, filtered) from the tracked .tdms file (both FP or TF tracker) '''
    t2xy = tracked2xy.Tracked_2_XY(trckd_file=trckd_file, 
                                   roi_num=roi_num, 
                                   plots_xy=plots_xy)
    xy2t = xy2torque.XY_2_Torque(t2xy.x, 
                                 t2xy.y,
                                 bead_diam_m=bead_diam_m,
                                 FPS=t2xy.FPS,
                                 umppx=umppx,
                                 rm_drift=rm_drift, 
                                 rm_drift_pts=rm_drift_pts, 
                                 rm_drift_mode=rm_drift_mode,
                                 rm_drift_plots=rm_drift_plots,
                                 stretch_xy=stretch_xy, 
                                 stretch_xy_plots=stretch_xy_plots,
                                 dist_beadsurf_wall=dist_beadsurf_wall,
                                 filter_name=filter_name, 
                                 filter_win=filter_win,
                                 plots=plots_torque)
    return xy2t
