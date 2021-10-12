#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Wed Dec 23 22:29:40 2020

@author: robertdenomme


jiggle trapezoid

Test
"""
import numpy as np
import math

import lagrangian as lagrangian
from lagrangian.dynamicalsystem import DynamicalSystem


qs_init = {'00': [0, 0],
           '10': [2, 0],
           '11': [3,2],
           '01': [1,2],
           }

qs_init = {i:np.array(qs_init[i]) for i in qs_init.keys()}

d = { i:{j:np.linalg.norm(np.array(qs_init[i]) - qs_init[j]) for j in qs_init.keys()}
     for i in qs_init.keys()}

def get_cos_theta(q, p_1, p_2):
    return np.dot(p_1 - q, p_2 - q) / (np.linalg.norm(p_1 - q) * np.linalg.norm(p_2 - q))


a = {'00': np.arccos(get_cos_theta(qs_init['00'], qs_init['10'], qs_init['01'])),
     '10': np.arccos(get_cos_theta(qs_init['10'], qs_init['00'], qs_init['11'])),
     '11': np.arccos(get_cos_theta(qs_init['11'], qs_init['10'], qs_init['01'])),
     '01': np.arccos(get_cos_theta(qs_init['01'], qs_init['00'], qs_init['11'])),
     }
     

masses = {i:1 for i in qs_init.keys()}

k_d = 100
k_a = 200

pi2 = math.pi/2

def gravity(*qs, g=9.8, masses=masses):
    u = 0
    for i in range(len(qs)):
        u += masses[list(qs_init.keys())[i]] * g * np.dot(qs[i] , [0,1])
    return u
        

def distance_constraint(q, p, *, rest_distance=1, spring_constant = 1):
    return spring_constant * (np.linalg.norm(q-p) - rest_distance)**2



def angle_constraint(q, p_1, p_2, *, rest_angle = math.pi/2, spring_constant = 1):
    cos_theta = get_cos_theta(q, p_1, p_2)
    return spring_constant * (cos_theta - math.cos(rest_angle))**2

        

dynamical_system = DynamicalSystem(qs_init,
                                   masses = masses,
                                   xlim = [-5,5],
                                   ylim = [-5, 5],
                                   wall_elasticity = .6,
                                   dt=1/60,
                                   integrator_code = 'ssprk3') #must have initial state to define the coordinates


dynamical_system.add_potential(distance_constraint,
                                ['00', '10'], 
                                rest_distance = d['00']['10'],
                                spring_constant = k_d
                                )

dynamical_system.add_potential(distance_constraint,
                                ['10', '11'], 
                                rest_distance = d['10']['11'],
                                spring_constant = k_d
                                )

dynamical_system.add_potential(distance_constraint,
                                ['11', '01'], 
                                rest_distance = d['11']['01'],
                                spring_constant = k_d
                                )

dynamical_system.add_potential(distance_constraint,
                                ['01', '00'], 
                                rest_distance = d['01']['00'],
                                spring_constant = k_d
                                )

'''Angle constraints'''


dynamical_system.add_potential(angle_constraint,
                                ['00', '10', '01'], 
                                rest_angle = a['00'],
                                spring_constant = k_a
                                )

dynamical_system.add_potential(angle_constraint,
                                ['10', '00', '11'], 
                                rest_angle = a['10'],
                                spring_constant = k_a
                                )

dynamical_system.add_potential(angle_constraint,
                                ['11', '10', '01'], 
                                rest_angle = a['11'],
                                spring_constant = k_a
                                )

dynamical_system.add_potential(angle_constraint,
                                ['01', '00', '11'], 
                                rest_angle = a['01'],
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
