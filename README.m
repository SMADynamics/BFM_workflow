# Bacterial Flagellar Motor analysis workflow 

This is the collection of the scripts we use to analyze BFM data from our bead assay experiments. 

Example:
    
    import tracked2torque
    
    # TF tracker: 
    # file points to the folder:
    file_to_open = '/home/francesco/ADYN/DATA/Amelie/240626/CL_240626_161019/'
    # or to 'trajectory.pt' file:
    file_to_open = '/home/francesco/CL_241202_144019/01/trajectory.pt'
    
    # FP tracker (Old): file should be the *_Trckd.tdms file name, eg:
    file_to_open = '/home/francesco/ADYN/DATA/Anais/201008/CL_201008_173841_Trckd.tdms'
    
    # run the full workflow, store all data in 'tt': 
    tt = tracked2torque.tracked_2_torque(
               trckd_file=file_to_open, 
               roi_num=1, 
               c0=1,                            # to crop the trace in [p0:p1]
               c1=None,     
               plots_xy=1,  
               bead_diam_m=1e-6,    
               dist_beadsurf_wall=10e-9,    
               umppx=0.1,                       # um to pixel
               correction_functions_order=['rm_outliers', 'rm_drift', 'stretch_xy'], 
               rm_outliers_plots=True,  
               rm_drift_mode='polyfit',         # ['spline' | 'linear' | 'polyfit']
               rm_drift_polyfit_deg=3,  
               rm_drift_wins=1001,              # win to define one interpolating point
               rm_drift_rmregions='',           # [(a1,b1), ..] or 'auto'. Remove bad regions of the trace (eg with pauses)
               rm_drift_rmregions_thr=0.8, 
               rm_drift_rmregions_extend=0.1,
               rm_drift_plots=True, 
               stretch_xy_plots=True,           # x,y ellipse fit and stretch to a circle
               filter_name='savgol',            # speed filter
               filter_win=101, 
               plots_torque=True)
    
    now in tt you have the following:
        'x_orig',
        'y_orig',
        'bead_diam_m',
        'FPS',
        'umppx',
        'correction_functions_order',
        'rm_outliers_findparam',
        'rm_outliers_win',
        'rm_outliers_plots',
        'rm_drift_mode',
        'rm_drift_wins',
        'rm_drift_qty',
        'rm_drift_qty_funct',
        'rm_drift_plots',
        'rm_drift_polyfit_deg',
        'rm_drift_rmregions',
        'rm_drift_rmregions_thr',
        'rm_drift_rmregions_extend',
        'stretch_xy_plots',
        'dist_beadsurf_wall_m',
        'filter_name',
        'filter_win',
        'plots',
        'plots_figname',
        'correction_functions',
        'x_corr',
        'y_corr',
        'angle_turns',
        'speed_Hz',
        'traj_radius_m',
        'drag_pNnms',
        'torque_pNnm',
        'torque_pNnm_filter',
        'speed_Hz_filter',
        'trckd_file',
        'c0',
        'c1',
        'workflow',
        'plots_all',
        'remove_outliers_funct',
        'remove_drift_funct',
        'stretch_xy_funct',

tracked2xy.py : extract from (FP or TF) tracked data (foldersand files) all the relevant information (x,y, metadata)

xy2torque.py : from xy data (coming from tracked2xy.py), calculate speed and torque.


Drift correction is controlled by the 'rm_drift_*' parameters:

    - rm_drift_mode : can be 'spline', 'linear' or 'polyfit'. Method to interpolate the points in the windows.
    - rm_drift_wins : numb. of time-windows to define one interpolating point for the drift correction function.
    - rm_drift_polyfit_deg : degree of the polynomial fit if rm_drift_mode is 'polyfit'.
    - rm_drift_qty='median' : quantity to calculate in each window, can be ['median'|'min'|'max'|'funct'].
    - rm_drift_qty_funct : function to use if rm_drift_qty is 'funct' (eg. np.mean))
    - rm_drift_rmregions : [(a1,b1), (a2,b2), ..] or 'auto'. Define bad regions (eg, pauses) of the trace to remove.
    - rm_drift_rmregions_thr : threshold if rm_drift_rmregions is 'auto'
    - rm_drift_rmregions_extend : extension of the bad regions if rm_drift_rmregions is 'auto' (eg, 0.1).
    - rm_drift_plots : make plots
