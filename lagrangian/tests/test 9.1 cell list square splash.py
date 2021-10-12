#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Wed Dec 23 22:29:40 2020

@author: robertdenomme


cell lists fluid molecules splash

"""
import numpy as np
import math
from itertools import product

import lagrangian as lagrangian
from lagrangian.dynamicalsystem import DynamicalSystem


N_rows = 8
N_columns = 12

N = N_rows * N_columns
r = 1


offset = np.array([-19.8,-19.1])

dt = 2**-9


velocity_init = np.array([0,0])
mass = 8

interatomic_distance = 1
cell_list_distance = 2
charge = 0



qs_init = dict()

for i,j in product(range(N_rows), range(N_columns)):
    qs_init[str(i)+'_'+str(j)] = np.array([r*i+ (r / 2) * j , r*j]) + offset


q_dots_init = {i:velocity_init for i in qs_init.keys()}

masses = {i:mass for i in qs_init.keys()}



def gravity(*qs, g=9.8, masses=masses):
    u = 0
    for i in range(len(qs)):
        u += masses[list(qs_init.keys())[i]] * g * np.dot(qs[i] , [0,1])
    return u
        

def van_der_waals(q, p, *, rest_distance=interatomic_distance, charge = charge):
    d = np.linalg.norm(q-p)
    eps = rest_distance
    vdw_potential = charge * d**-1 + eps*((d/eps)**-12 - 2 * (d/eps)**-6)
    return vdw_potential



     

dynamical_system = DynamicalSystem(qs_init, 
                                   q_dots_init,
                                   masses = masses,
                                   xlim = [-20,20],
                                   ylim = [-20, 20],
                                   cell_list_dx = cell_list_distance,
                                   cell_list_dy = cell_list_distance,
                                   wall_elasticity = .6,
                                   dt=dt,
                                   integrator_code = 'semi-implicit euler') #must have initial state to define the coordinates

for q in dynamical_system.initial_state.qs:
    dynamical_system.add_rendered_path([q] )



dynamical_system.add_cell_list_pair_potential(van_der_waals)


'''gravity'''

dynamical_system.add_potential(gravity, #the potential function
                                [name for name in qs_init.keys()], #the argument coordinate names
                                g=9.8, masses=masses #the keyword arguments to pass
                                )


time = 5

dynamical_system.run_dynamics(time)


dynamical_system.display()
