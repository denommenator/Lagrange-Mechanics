#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Sun Dec 27 17:15:39 2020

@author: robertdenomme
"""


from lagrangian.state import State, Coordinates



def resolve_wall_collisions(state, 
                            xlim = (-10, 10), 
                            ylim = (-10, 10),
                            wall_elasticity = 1):
    new_state = state.copy()
    
    eps = .0001
    e = wall_elasticity
    
    qs = new_state.qs
    q_dots = new_state.q_dots
    
    for i, q in qs.items():
        if q[0] < xlim[0]:
            qs[i][0] = xlim[0]+eps
            if q_dots[i][0] < 0:
                q_dots[i][0] = - e * q_dots[i][0]
        
        elif q[0] > xlim[1]:
            qs[i][0] = xlim[1]-eps
            if q_dots[i][0] > 0:
                q_dots[i][0] = - e * q_dots[i][0]
        
        
        if q[1] < ylim[0]:
            qs[i][1] = ylim[0]+eps
            if q_dots[i][1] < 0:
                q_dots[i][1] = - e * q_dots[i][1]
        
        elif q[1] > ylim[1]:
            qs[i][1] = ylim[1]-eps
            if q_dots[i][1] > 0:
                q_dots[i][1] = - e * q_dots[i][1]
    
    return State(qs, q_dots)