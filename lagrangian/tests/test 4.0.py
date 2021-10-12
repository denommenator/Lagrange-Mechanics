#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Wed Dec 23 22:29:40 2020

@author: robertdenomme

elastic band dropping and bouncing off the floor


Test
"""
import numpy as np

import lagrangian as lagrangian
from lagrangian.dynamicalsystem import DynamicalSystem


N = 30
dx = 1/N
qs_init = {'q'+str(i): [i*dx, 0] for i in range(N)}


d = [ np.linalg.norm(np.array(qs_init['q'+str(i+1)])-qs_init['q'+str(i)]) for i in range(N-1)]
# d_12 = np.linalg.norm(np.array(qs_init['q1'])-qs_init['q2'])
# d_23 = np.linalg.norm(np.array(qs_init['q2'])-qs_init['q3'])
# d_34 = np.linalg.norm(np.array(qs_init['q3'])-qs_init['q4'])
# d_45 = np.linalg.norm(np.array(qs_init['q4'])-qs_init['q5'])
                               
masses = {'q'+str(i):10*1/N for i in range(N)}

def gravity(*qs, g=9.8, masses=masses):
    u = 0
    for i in range(len(qs)):
        u += masses['q'+str(i)] * g * np.dot(qs[i] , [0,1])
    return u
        

def distance_constraint(q, p, *, rest_distance=0):
    return 100 * (np.linalg.norm(q-p) - rest_distance)**2

dynamical_system = DynamicalSystem(qs_init,
                                   masses = masses,
                                   xlim = [0,1],
                                   ylim = [-2.5, 1],
                                   wall_elasticity = 2,
                                   dt=1/60,
                                   integrator_code = 'ssprk3') #must have initial state to define the coordinates


dynamical_system.add_potential(distance_constraint,
                                ['q0'],
                                p=[0,0] #pass p as a keyword argument to fix distance from q1 to origin
                                )

dynamical_system.add_potential(distance_constraint,
                                ['q'+str(N-1)],
                                p=qs_init['q'+str(N-1)] #pass p as a keyword argument to fix distance from q1 to origin
                                )

for i in range(N-1):
    dynamical_system.add_potential(distance_constraint,
                                    ['q'+str(i),'q'+str(i+1)],
                                    rest_distance = d[i]
                                    )


dynamical_system.add_potential(gravity, #the potential function
                               ['q'+str(i) for i in range(N)], #the argument coordinate names
                               g=9.8, masses=masses #the keyword arguments to pass
                               )


time = 5

dynamical_system.run_dynamics(time)


dynamical_system.display()
