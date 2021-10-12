#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Wed Dec 23 22:29:40 2020

@author: robertdenomme


Test
"""
import numpy as np
import math
import itertools

import lagrangian as lagrangian
from lagrangian.dynamicalsystem import DynamicalSystem



def stable_orbit_velocity(M, d, G):
        return math.sqrt(G*M*d**-1)
    



N = 4
0

G = 16

dt = 2**-5


d_earth = 372
d_moon = 1

qs_init = {
    'sun': [0,0],
    'earth': [d_earth, 0],
    'moon': [d_earth + d_moon, 0]
    }

masses = {
    'sun': 330_000,
    'earth': 1,
    'moon': .25
    }


v_earth = stable_orbit_velocity(masses['sun'],d_earth, G=G)
v_moon = stable_orbit_velocity(masses['earth'],d_moon, G=G)

q_dots_init = {
    'sun': [0,0],
    'earth': [0,v_earth],
    'moon': [0,v_earth + 1 * v_moon]
    }







def gravity(q, p, mass1, mass2, G=G):
    d = np.linalg.norm(q-p)
    
    
    return - mass1 * mass2 * G / d
        



     

dynamical_system = DynamicalSystem(qs_init, 
                                   q_dots_init,
                                   masses = masses,
                                   xlim = [-400,400],
                                   ylim = [-400, 400],
                                   # wall_collisions = False,
                                   dt=dt,
                                   integrator_code = 'ssprk3') #must have initial state to define the coordinates

for i in qs_init.keys():
    dynamical_system.add_rendered_path([i] )



for i, j in itertools.combinations(qs_init.keys(), 2):
    dynamical_system.add_potential(gravity,
                                [i, j], 
                                mass1 = masses[i],
                                mass2 = masses[j]
                                )




time = 10

dynamical_system.run_dynamics(time)


dynamical_system.display()
