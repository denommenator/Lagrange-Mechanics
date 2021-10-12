#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Mon Dec 21 17:47:39 2020

@author: robertdenomme

The main job of the phase space is the function that computes forces.
The phase space will also be responsible for running the dynamics, i.e.
passing the forces to the integrator and adding the next state into a 
trajectories object.


"""

import lagrangian.integrators as integrators

import numpy as np
from itertools import product
    
integrator_dict = integrators.integrator_dict

class PhaseSpace:
    def __init__(self, particle_system,
                 *,
                 xlim = (-10,10), 
                 ylim = (-10,10), 
                 drag = 0,
                 g = 9.8,
                 dt = 1/60,
                 integrator_code = 'ssprk3'):
        
        self.xlim = xlim
        self.ylim = ylim
        self.drag = drag
        self.g = g
        self.dt = dt
        self.integrator = integrator_dict[integrator_code]
        
        
    
    