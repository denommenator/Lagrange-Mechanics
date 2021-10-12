#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Thu Dec 24 23:37:20 2020

@author: robertdenomme
"""


import matplotlib.pyplot as plt
import matplotlib.animation as animation

from matplotlib.path import Path
import matplotlib.patches as patches


from matplotlib.path import Path
import matplotlib.patches as patches

import numpy as np

class MatPlotRenderer:
    def __init__(self, trajectory_data):
        self.trajectory_data = trajectory_data
        self.dynamical_system = trajectory_data.dynamical_system
        
        
        self.xlim = self.dynamical_system.xlim
        self.ylim = self.dynamical_system.ylim
        
        
        self.figsize = (10,10)
        
        self.__setup__()
        
    def __setup__(self):
        
        self.fig = plt.figure(figsize=self.figsize)
        self.ax = self.fig.add_subplot(111, 
                                       autoscale_on=False, 
                                       xlim=self.xlim, 
                                       ylim=self.ylim)
        self.ax.set_aspect('equal')
        self.ax.grid()


        self.time_template = 'time = %.1fs'
        self.time_text = self.ax.text(0.05, 0.9, '', transform=self.ax.transAxes)
        
        
        
        
        
        self.trajectory_data.process_data_for_rendering()
        
        
        
        self.verts = np.array(self.trajectory_data.vertices_at_time[0])
        self.codes = self.dynamical_system.rendered_path_codes
        
        self.paths = Path(self.verts, self.codes)

        self.patch = patches.PathPatch(self.paths, 
                                  lw=2,
                                  animated = True)
        self.ax.add_patch(self.patch)

                


    def __animate__(self, t):
        
        self.verts[:] = self.trajectory_data.vertices_at_time[t]
        self.time_text.set_text(self.time_template % (t * self.trajectory_data.frame_dt))
        
        x,y = zip(*self.paths.vertices)
        dots, = self.ax.plot(x, y, 'ro', markersize=8)
        
        return  dots, self.patch, self.time_text,


    

    def display(self):
        
        
        self.ani = animation.FuncAnimation(self.fig, 
                                          self.__animate__, 
                                          self.trajectory_data.N_frames, 
                                          interval=self.trajectory_data.frame_dt * 1000, 
                                          #interval = 1000,
                                          blit=True)
        
        plt.show()


