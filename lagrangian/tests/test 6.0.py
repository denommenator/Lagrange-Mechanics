#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Wed Dec 23 22:29:40 2020

@author: robertdenomme


Squishy ball
"""
import numpy as np
import math

import lagrangian as lagrangian
from lagrangian.dynamicalsystem import DynamicalSystem


N = 20
r = 2
qs_init = {str(i): [r*math.cos(i * 2*math.pi / N), r*math.sin(i * 2*math.pi / N)]
           for i in range(N)}

offset = np.array([0,-.0001])
qs_init = {i:np.array(qs_init[i]) + offset for i in qs_init.keys()}

velocity_init = np.array([4,0])
q_dots_init = {i:velocity_init for i in qs_init.keys()}

N = len(qs_init)

d = { str(i):np.linalg.norm(np.array(qs_init[str((i+1) % N)]) - qs_init[str(i % N)]) \
     for i in range(N)}



    
# def get_cos_theta(q, p_1, p_2):
#     return np.dot(p_1 - q, p_2 - q) / (np.linalg.norm(p_1 - q) * np.linalg.norm(p_2 - q))


# def get_orientation_sign(q, p_1, p_2):
#     '''Need more than just the cosine of the interior angle between the two.
#     '''
#     v = p_1 - q
#     w = p_2 - q
    
#     a, b = v[0], v[1]
#     c, d = w[0], w[1]
    
#     orientation_sign = np.sign(a * d - b * c)
    
#     return orientation_sign


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


# for i in range(N):
#     print(i, ': ',
#         get_theta(qs_init[str(i % N)], qs_init[str((i-1) % N)], qs_init[str((i+1) % N)]),
#         '\n'
#         )

a = {str(i): get_theta(qs_init[str(i % N)], qs_init[str((i-1) % N)], qs_init[str((i+1) % N)])
    for i in range(N)
     }
     

masses = {i:1*1/N for i in qs_init.keys()}

k_d = 4 * N
#distance spring coefficient

k_a = 2**-1 * N
#angular spring coefficient


def gravity(*qs, g=9.8, masses=masses):
    u = 0
    for i in range(len(qs)):
        u += masses[list(qs_init.keys())[i]] * g * np.dot(qs[i] , [0,1])
    return u
        

def distance_constraint(q, p, *, rest_distance=1, spring_constant = 1):
    return spring_constant * (np.linalg.norm(q-p) - rest_distance)**2



def angle_constraint(q, p_1, p_2, *, rest_angle = math.pi/2, spring_constant = 1):
    theta = get_theta(q, p_1, p_2)
    return spring_constant * (theta - rest_angle)**2

        

dynamical_system = DynamicalSystem(qs_init, 
                                   q_dots_init,
                                   masses = masses,
                                   xlim = [-5,5],
                                   ylim = [-5, 5],
                                   wall_elasticity = 1,
                                   dt=1/2**8,
                                   integrator_code = 'verlet') #must have initial state to define the coordinates

dynamical_system.add_rendered_path(
                                    [#one path
                                    str(i) for i in range(N)
                                    ] + [str(0)]
                                        )


for i in range(N):
    dynamical_system.add_potential(distance_constraint,
                                [str(i), str((i+1) % N)], 
                                rest_distance = d[str(i)],
                                spring_constant = k_d
                                )

'''Angle constraints'''

for i in range(N):
    dynamical_system.add_potential(angle_constraint,
                                [str(i), str((i-1) % N), str((i+1) % N)], 
                                rest_angle = a[str(i)],
                                spring_constant = k_a
                                )



'''gravity'''

dynamical_system.add_potential(gravity, #the potential function
                               [name for name in qs_init.keys()], #the argument coordinate names
                               g=9.8, masses=masses #the keyword arguments to pass
                               )


time = 10

dynamical_system.run_dynamics(time)


dynamical_system.display()
