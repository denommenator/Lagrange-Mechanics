#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Wed Dec 23 22:29:40 2020

@author: robertdenomme


Test
"""
import numpy as np

import lagrangian as lagrangian
from lagrangian.dynamicalsystem import DynamicalSystem

qs_init = {'q1':[0,0],
           'q2':[3,4],
           'q3':[5,6],
           }


d_12 = np.linalg.norm(np.array(qs_init['q1'])-qs_init['q2'])
d_23 = np.linalg.norm(np.array(qs_init['q2'])-qs_init['q3'])

                               

def gravity(*qs, g=9.8, masses=[1,1,1]):
    u = 0
    for i, q_i in enumerate(qs):
        u += masses[i] * g * np.dot(q_i , [0,1])
    return u
        

def distance_constraint(q, p, *, rest_distance=0):
    return 100 * (np.linalg.norm(q-p) - rest_distance)**2

dynamical_system = DynamicalSystem(qs_init) #must have initial state to define the coordinates




dynamical_system.add_potential(distance_constraint,
                               ['q1','q2'],
                               rest_distance = d_12,
                               )

dynamical_system.add_potential(distance_constraint,
                               ['q2','q3'],
                               rest_distance = d_23
                               )

dynamical_system.add_potential(distance_constraint,
                               ['q1'],
                               p=[0,0] #pass p as a keyword argument to fix distance from q1 to origin
                               )

dynamical_system.add_potential(gravity, #the potential function
                               ['q1', 'q2', 'q3'], #the argument coordinate names
                               g=9.8, masses=[1,1,1] #the keyword arguments to pass
                               )


time = 40

dynamical_system.run_dynamics(time)


dynamical_system.display()
