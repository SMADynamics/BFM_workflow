import numpy as np
import matplotlib.pyplot as plt
import torch
import json
import re
import os

import openTDMS





class Tracked_2_XY():
    
    def __init__(self, trckd_file='', 
                 roi_num=0, 
                 c0=0, c1=None, 
                 plots_xy=False):
        ''' get x,y of the bead fropm tracked file (FP or TF tracker)'''
        self.trckd_file = trckd_file
        self.roi_num = roi_num
        self.c0 = c0
        self.c1 = c1
        print(f'Tracked_2_XY.__init__(): trckd_file: {trckd_file}; roi_num: {roi_num}; c0: {c0}; c1: {c1}')
        self.workflow(plots_xy=plots_xy)
        print(f'Tracked_2_XY.__init__(): Done.')


    def workflow(self, plots_xy=False):
        self.get_trck()
        self.get_roi(roi_num=self.roi_num, c0=self.c0, c1=self.c1)
        if plots_xy:
            self.plots_xy()


    def get_trck(self):
        '''Get tracked data from trckd_file, from
            FP tracker ('/CL*_Trckd.tdms') or
            TF tracker ('/trajectory.pt') '''
        # FP tracker :
        if self.trckd_file.endswith('_Trckd.tdms'):
            print(f'Tracked_2_XY.get_trck(): FP-tracked file: {self.trckd_file}')
            self.trckd = openTDMS.openTdmsFile(self.trckd_file)
            print(f'Tracked_2_XY.get_trck(): keys found: {self.trckd.keys()}')
            nrois = re.search(r'Number of ROI : (\d)', self.trckd['/CL-config/#X'][0]).group(1)
            print(f'Tracked_2_XY.get_trck(): Num. or ROIs: {nrois}')
        # TF tracker:
        else:
            # if .tdms file provided:
            if self.trckd_file.endswith('.tdms'):
                directory     = os.path.split(self.trckd_file)[0]
                file          = os.path.split(self.trckd_file)[1].removesuffix('.tdms')
                traj_file     = os.path.join(directory, file, str(self.roi_num), 'trajectory.pt')
                metadata_file = os.path.join(directory, file, 'metadata.json')
            # if trajectory.pt file provided:
            elif self.trckd_file.endswith('trajectory.pt'):
                m = re.search(r'^(.*?/CL_[^/]+/)', self.trckd_file)
                directory = m.group(1) if m else None
                traj_file = self.trckd_file
                metadata_file = os.path.join(directory, 'metadata.json')
            print(f'Tracked_2_XY.get_trck(): TF-tracker traj_file: {traj_file}')
            print(f'Tracked_2_XY.get_trck(): TF-tracker directory: {directory}')
            print(f'Tracked_2_XY.get_trck(): TF-tracker metadata_file: {metadata_file}')
            self.trckd = torch.load(traj_file, weights_only=True)
            self.metadata = json.load(open(metadata_file))
            print(f'Tracked_2_XY.get_trck(): TF-tracker file loaded.')
            #print(f'Tracked_2_XY.get_trck(): TF-tracker metadata: {self.metadata}')


    def get_roi(self, c0=0, c1=-1, roi_num=0):
        ''' get info (x,y,FPS) from ROI number 'roi_num' from tracked file
            (guessing FP or TF tracker)
        '''
        # FP tracker file:
        if self.trckd_file.endswith('_Trckd.tdms'):
            if not hasattr(self, 'trckd'):
                self.get_trck()
            # get x,y:
            x_key = f'/ROI{roi_num}_Trk/X{roi_num}'
            y_key = f'/ROI{roi_num}_Trk/Y{roi_num}'
            self.x = np.float64(self.trckd[x_key][c0:c1])
            self.y = np.float64(self.trckd[y_key][c0:c1])
            print(f'Tracked_2_XY.get_roi(): roi_num:{roi_num}; found x,y of {len(self.x)} pts; type: {self.x.dtype}; interval [c0:c1]: [{c0}:{c1}]')
            # get FPS:
            fps_idx_start = self.trckd['/CL-config/#X'][0].find('Frame Rate : ') + len('Frame Rate : ')
            fps_idx_end = self.trckd['/CL-config/#X'][0].find('\r', fps_idx_start)
            self.FPS = float(self.trckd['/CL-config/#X'][0][fps_idx_start:fps_idx_end])
            print(f'Tracked_2_XY.get_roi(): found FPS {self.FPS}')
        # TF tracker files:
        else:
            #self.trckd is a tensor in a dictionary with keys cx, cy and z. 
            # To float64 to avoid problems unwrapping the angle:
            self.x = np.float64(self.trckd['cy'].numpy()[c0:c1])
            self.y = np.float64(self.trckd['cx'].numpy()[c0:c1])
            self.FPS = float(self.metadata['Frame Rate'])
            print(f'Tracked_2_XY.get_roi(): populated with x(t), y(t) (len: {len(self.x)} type: {self.x.dtype}), FPS={self.FPS}')


    def plots_xy(self):
        fig = plt.figure('Tracked_2_XY.plots_xy', clear=True)
        ax1 = fig.add_subplot(311)
        ax2 = fig.add_subplot(312)
        ax3 = fig.add_subplot(313)
        ax1.plot(self.x, self.y, ',', alpha=0.2)
        ax1.axis('equal')
        ax1.set_xlabel('x (px)')
        ax1.set_ylabel('y (px)')
        ax2.plot(self.x, ',', alpha=0.2)
        ax2.set_ylabel('x (px)')
        ax2.set_xlabel('idx')
        ax3.plot(self.y, ',', alpha=0.2)
        ax3.set_ylabel('y (px)')
        ax3.set_xlabel('idx')
        # Create a secondary x-axis for time
        def index_to_time(x):
            return x / self.FPS  # Convert index to time
        def time_to_index(x):
            return x * self.FPS  # Convert time back to index
        sec_ax2 = ax2.secondary_xaxis('top', functions=(index_to_time, time_to_index))
        sec_ax3 = ax3.secondary_xaxis('top', functions=(index_to_time, time_to_index))
        sec_ax2.set_xlabel("Time (seconds)")
        sec_ax3.set_xlabel("Time (seconds)")
        fig.tight_layout()






