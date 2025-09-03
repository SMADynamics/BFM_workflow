# calls to get torque (and speed, filtered) from the tracked .tdms file (both FP or TF tracker)

import tracked2xy
import xy2torque



def tracked_2_torque(trckd_file='', 
                     roi_num=0,
                     c0=0,
                     c1=None,
                     plots_xy=False,
                     bead_diam_m=1e-6, 
                     dist_beadsurf_wall_m=100e-6,
                     umppx=0.1,
                     correction_functions_order=['rm_outliers', 'rm_drift', 'stretch_xy'],
                     rm_outliers_findparam=3,
                     rm_outliers_win=5,
                     rm_outliers_plots=False,
                     rm_drift_mode='linear',
                     rm_drift_pts=100, 
                     rm_drift_plots=False,
                     stretch_xy_plots=False,
                     filter_name='median', 
                     filter_win=101,
                     plots_torque=False,
                     store_corr=False ):
    ''' 
        calls to get torque (and speed, filtered) from the tracked .tdms file (both FP or TF tracker) 
        
        Example:
          
            > file = '/home/francesco/ADYN/DATA/Amelie/240626/CL_240626_161019/'
            > tracked2torque.tracked_2_torque(trckd_file=file, roi_num=1, c0=1, c1=None, plots_xy=1, 
                        bead_diam_m=1e-6, dist_beadsurf_wall_m=10e-9, 
                        umppx=0.1, 
                        correction_functions_order=['rm_outliers', 'rm_drift', 'stretch_xy'], 
                        rm_outliers_plots=True, rm_drift_mode='spline', rm_drift_plots=True, 
                        stretch_xy_plots=True, 
                        filter_name='savgol', filter_win=101, 
                        plots_torque=True)
    '''
    tr2xy = tracked2xy.Tracked_2_XY(trckd_file=trckd_file, 
                                   roi_num=roi_num, 
                                   c0=c0, 
                                   c1=c1,
                                   plots_xy=plots_xy)
    xy2tq = xy2torque.XY_2_Torque(tr2xy.x, 
                                 tr2xy.y,
                                 bead_diam_m=bead_diam_m,
                                 FPS=tr2xy.FPS,
                                 umppx=umppx,
                                 correction_functions_order=correction_functions_order,
                                 rm_outliers_findparam=rm_outliers_findparam,
                                 rm_outliers_win=rm_outliers_win,
                                 rm_outliers_plots=rm_outliers_plots,
                                 rm_drift_pts=rm_drift_pts, 
                                 rm_drift_mode=rm_drift_mode,
                                 rm_drift_plots=rm_drift_plots,
                                 stretch_xy_plots=stretch_xy_plots,
                                 dist_beadsurf_wall_m=dist_beadsurf_wall_m,
                                 filter_name=filter_name, 
                                 filter_win=filter_win,
                                 plots=plots_torque,
                                 store_corr=store_corr)
    return xy2tq
