#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Wed Dec 23 22:29:40 2020

@author: robertdenomme


fluid molecules splash

"""
import numpy as np
import math
import itertools

import lagrangian as lagrangian
from lagrangian.dynamicalsystem import DynamicalSystem


N = 4
r = .5


offset = np.array([0,-.0001])

dt = 2**-9

velocity_init = np.array([4,0])
mass = 2**-2

interatomic_distance = 1
charge = 0


phi = (1 + math.sqrt(5)) / 2
theta_0 = 2 * math.pi / phi




qs_init = dict()

for n in range(0, N):
    position = np.array([n**.5 * math.cos(n*theta_0), n**.5 * math.sin(n*theta_0)])
    normalized_position = r * position
    qs_init[str(n)] = normalized_position + offset


qs_init['0'] = np.array([0,0])
qs_init['1'] = np.array([3,0])


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

def distance_constraint(q, p, *, rest_distance=1, spring_constant = 1):
    return spring_constant * (np.linalg.norm(q-p) - rest_distance)**2



     

dynamical_system = DynamicalSystem(qs_init, 
                                   q_dots_init,
                                   masses = masses,
                                   xlim = [-5,5],
                                   ylim = [-5, 5],
                                   wall_elasticity = 1,
                                   dt=dt,
                                   integrator_code = 'verlet') #must have initial state to define the coordinates

for i in range(N):
    dynamical_system.add_rendered_path([str(i)] )



for i, j in itertools.combinations(range(N), 2):
    dynamical_system.add_potential(van_der_waals,
                                [str(i), str(j)], 
                                )


'''gravity'''

dynamical_system.add_potential(gravity, #the potential function
                                [name for name in qs_init.keys()], #the argument coordinate names
                                g=9.8, masses=masses #the keyword arguments to pass
                                )


time = 10

dynamical_system.run_dynamics(time)


dynamical_system.display()
