# Bacterial Flagellar Motor analysis workflow

This is the collection of the scripts we use to analyze BFM data from our bead assay experiments. 

Example:
    # file should be the folder for TF tracker, eg:
    file_to_open = '/home/francesco/ADYN/DATA/Amelie/240626/CL_240626_161019/'
    # or the *_Trckd.tdms file name for FP tracker, eg:
    file_to_open = '/home/francesco/ADYN/DATA/Anais/201008/CL_201008_173841_Trckd.tdms'
    # run the full workflow: 
    tracked2torque.tracked_2_torque(
               trckd_file=file_to_open, 
               roi_num=1, 
               c0=1, c1=None, 
               plots_xy=1, 
               bead_diam_m=1e-6, 
               dist_beadsurf_wall=10e-9, 
               umppx=0.1, 
               correction_functions_order=['rm_outliers', 'rm_drift', 'stretch_xy'], 
               rm_outliers_plots=True, 
               rm_drift_mode='spline', 
               rm_drift_plots=True, 
               stretch_xy_plots=True, 
               filter_name='savgol', 
               filter_win=101, 
               plots_torque=True)


tracked2xy.py : extract from (FP or TF) tracked data (foldersand files) all the relevant information (x,y, metadata)

xy2torque.py : from xy data (coming from tracked2xy.py), calculate speed and torque.

