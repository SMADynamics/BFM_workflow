import numpy as np
import matplotlib.pyplot as plt
import torch
import json
import re
import os

import openTDMS





class Tracked_2_XY():
    
    def __init__(self, trckd_file='', roi_num=0, c0=0, c1=None, plots_xy=False):
        ''' get x,y of the bead fropm tracked file (FP or TF tracker)'''
        self.trckd_file = trckd_file
        self.roi_num = roi_num
        self.c0 = c0
        self.c1 = c1
        self.workflow(plots_xy=plots_xy)



    def workflow(self, plots_xy=False):
        self.get_trck()
        self.get_roi(roi_num=self.roi_num, c0=self.c0, c1=self.c1)
        if plots_xy:
            self.plots_xy()



    def get_trck(self):
        '''Get tracked data from trckd_file, from
            FP tracker ('_Trckd.tdms') or
            TF tracker ('trajectory.pt') '''
        # FP tracker:
        if self.trckd_file.endswith('_Trckd.tdms'):
            print(f'get_trck(): FP-tracked file: {self.trckd_file}')
            self.trckd = openTDMS.openTdmsFile(self.trckd_file)
            print(f'get_trck(): keys found: {self.trckd.keys()}')
            nrois = re.search(r'Number of ROI : (\d)', self.trckd['/CL-config/#X'][0]).group(1)
            print(f'get_trck(): Num. or ROIs: {nrois}')
        # TF tracker:
        else:
            directory = os.path.split(self.trckd_file)[0]
            file = os.path.split(self.trckd_file)[1].strip('.tdms')
            self.trckd = torch.load(directory + '/' + file + '/' + str(self.roi_num) + '/trajectory.pt')
            metadata_file = open(directory + '/' + file + '/' + 'metadata.json')
            self.metadata = json.load(metadata_file)
            print(f'get_trck(): TF-tracker file loaded.')
            #print(f'get_trck(): TF-tracker metadata: {self.metadata}')



    def get_roi(self, c0=0, c1=-1, roi_num=1):
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
            self.x = self.trckd[x_key][c0:c1]
            self.y = self.trckd[y_key][c0:c1]
            print(f'get_roi(): roi_num:{roi_num}; found x,y of {len(self.x)} pts; interval [c0:c1]: [{c0}:{c1}]')
            # get FPS:
            fps_idx_start = self.trckd['/CL-config/#X'][0].find('Frame Rate : ') + len('Frame Rate : ')
            fps_idx_end = self.trckd['/CL-config/#X'][0].find('\r', fps_idx_start)
            self.FPS = float(self.trckd['/CL-config/#X'][0][fps_idx_start:fps_idx_end])
            print(f'get_roi(): found FPS {self.FPS}')
        # TF tracker files:
        else:
            #self.trckd is a tensor in a dictionary with keys cx, cy and z.
            self.x = self.trckd['cy'].numpy()[c0:c1]
            self.y = self.trckd['cx'].numpy()[c0:c1]
            self.FPS = float(self.metadata['Frame Rate'])
            print(f'get_roi(): popultaed with x(t), y(t) (len: {len(self.x)}), and FPS={self.FPS}')


    def plots_xy(self):
        fig = plt.figure('plots_xy', clear=True)
        ax1 = fig.add_subplot(321)
        ax2 = fig.add_subplot(312)
        ax3 = fig.add_subplot(313)
        ax1.plot(self.x, self.y, ',', alpha=0.2)
        ax1.set_xlabel('x (px)')
        ax1.set_ylabel('y (px)')
        ax2.plot(self.x, ',', alpha=0.4)
        ax2.set_ylabel('x (px)')
        ax2.set_xlabel('idx')
        ax3.plot(self.y, ',', alpha=0.4)
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






