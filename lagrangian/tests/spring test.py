#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Wed Dec 30 11:34:49 2020

@author: robertdenomme

Linear and angular spring stability test suite. 
Used for testing stability of integrators on the simply harmonic osciallator 
at different energy levels and with different masses
"""

import numpy as np
import math

import lagrangian as lagrangian
from lagrangian.dynamicalsystem import DynamicalSystem


N_springs = 6


longest_spring_length = 4

linear_spring_constant = 2 * 5
linear_spring_distance = 1
mass = 1/5



l = longest_spring_length
qs_init = dict()
for i in range(N_springs):
    qs_init['left_spr_'+str(i)] = [-l * 2**-1 * 2**(-i), i * .3] # spring i has length 2 * 2**-i
    qs_init['right_spr_'+str(i)] = [l * 2**-1 * 2**(-i), i * .3]
           

offset = np.array([0,-.0001])
qs_init = {i:np.array(qs_init[i]) + offset for i in qs_init.keys()}

masses = {i:mass for i in qs_init.keys()}





def get_theta(q, p_1, p_2):
    v = p_1 - q
    w = p_2 - q
    
    a, b = v[0], v[1]
    norm_v = np.linalg.norm(v)
    
    c, d = w[0], w[1]
    norm_w = np.linalg.norm(w)
    
    det = a * d - b * c
    sin_theta = det / (norm_v * norm_w)
    
    
    
    return np.arcsin(sin_theta)









def distance_constraint(q, p, *, rest_distance=1, spring_constant = 1):
    return spring_constant * (np.linalg.norm(q-p) - rest_distance)**2



def angle_constraint(q, p_1, p_2, *, rest_angle = math.pi/2, spring_constant = 1):
    theta = get_theta(q, p_1, p_2)
    return spring_constant * (theta - rest_angle)**2

        

dynamical_system = DynamicalSystem(qs_init, 
                                   
                                   masses = masses,
                                   xlim = [-10,10],
                                   ylim = [-1, N_springs * .3 +1],
                                   wall_elasticity = 1,
                                   dt=1/60,
                                   integrator_code = 'verlet') #must have initial state to define the coordinates

for s in range(N_springs):
    dynamical_system.add_rendered_path(
                                    [
                                    'left_spr_'+str(s),
                                    'right_spr_'+str(s)
                                    ] 
                                        )


for s in range(N_springs):
    dynamical_system.add_potential(distance_constraint,
                                ['left_spr_'+str(s), 'right_spr_'+str(s)], 
                                rest_distance = linear_spring_distance,
                                spring_constant = linear_spring_constant
                                )


time = 10

dynamical_system.run_dynamics(time)


dynamical_system.display()
