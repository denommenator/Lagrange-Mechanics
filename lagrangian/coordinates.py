#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Mon Dec 21 16:40:36 2020

@author: robertdenomme

coordinates and state are objects that should hold only the data for a 
system that changes from time step to time step. The dynamics engine 
(phase-space) will need this data, along with the data that stays constant 
 from time step to time step (a dynamical system)
"""


import numpy as np


class Coordinate: #???
    def __init__(self, coords):
        self.__coords = np.array(coords)
        
        
class State:
    def __init__(self):
        
        self.qs = dict()
        self.q_dots = dict()
        
    
    def add_coordinate(self, name = None, value = None):
        if name == None:
            name = len(self.coordinates.keys())
        
        self.coordinates[name] = Coordinate(value)
    
    def set_coordinate(self, name, value):
        self.coordinates[name] = Coordinate(value)
        
        
    