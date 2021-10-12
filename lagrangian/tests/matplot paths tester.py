#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Tue Dec 29 22:55:46 2020

@author: robertdenomme
"""

#matplot tester


import numpy as np


import matplotlib.pyplot as plt
import matplotlib.animation as animation

from matplotlib.path import Path
import matplotlib.patches as patches



fig = plt.figure(figsize=(10,10))
ax = fig.add_subplot(111, 
                          autoscale_on=False, 
                          xlim = (0,10),
                          ylim = (0,10)
                          )
ax.set_aspect('equal')
ax.grid()


time_template = 'time = %.1fs'
time_text = ax.text(0.05, 0.9, '', transform=ax.transAxes)
        

verts = np.array([(0,0),
                (1,2),
                (3,4),
                (4,4),
                (5,3),
                (6,2),
                (0,0) #ignored
                ], dtype = float)

codes = [Path.MOVETO] + [Path.LINETO]*(len(verts)-2) + [Path.CLOSEPOLY]


path = Path(verts, codes)

patch = patches.PathPatch(path, 
                          facecolor = 'orange', 
                          lw=2,
                          animated = True)
ax.add_patch(patch)


#artists = ax.plot([0,1,2],[1,2,3], 'ro')

def animate(t):
    
    verts[:] = verts + np.array([.1,.1])
    time_text.set_text(time_template % (t * 1/20))
    
    return  patch, time_text,



ani = animation.FuncAnimation(fig, 
                                          animate, 
                                          20 * 5, 
                                          interval= 1000 / 20, 
                                          #interval = 1000,
                                          blit=True)
        
plt.show()
        
        
        