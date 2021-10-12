#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Thu Dec 24 23:49:32 2020

@author: robertdenomme
"""


import numpy as np

from lagrangian.state import State, Coordinates


def verlet_next(state, dt = None):
    if dt is None:
        dt = state.dynamical_system.dt
        
    current_state = state.qs
    previous_state = state.dynamical_system.trajectory_data[-2]
    
    q_dotdots = state.get_acceleration()
    
    qs_next = 2 * state.qs - previous_state.qs + dt**2  * q_dotdots 
    
    q_dots_next = state.q_dots #(we aren't updating these because they are not needed. If they are needed, use velocity verlet
    
    state_next = State(qs_next, q_dots_next)
    return state_next

    
def deprecated_verlet_next(state, phase_space, dt=None):
    '''We're going to assume q_n, and q_dots_n are related as follows
    
    q_n+1 = 2q_n + q_dot_n dt + 1/2 q_dotdot dt^2
    
    q_dot_n = (q_n+1 - 2 q_n)/dt - 1/2 q_dotdot dt
    
    
    giving the relation 
    
    q_n+1 = 2q_n - q_n-1 + q_dotdot dt^2
    
    '''
    if dt == None:
        dt = phase_space.dt
    
    forces = phase_space.forces(state)
    q_dotdots = phase_space.masses **-1 * forces
    
    
    if len(phase_space.trajectories) == 1:
        qs_next = 2 * state.qs + state.q_dots * dt + 1/2 * q_dotdots * dt**2
    else: 
        prev_state = phase_space.trajectories[-2]
        qs_next = 2*state.qs - prev_state.qs + q_dotdots * dt**2
    
    q_dots_next = (qs_next- 2*state.qs)/dt - 1/2 * q_dotdots * dt
    
    state_next = State(qs_next, q_dots_next)
    return state_next


def forward_euler_next(state, dt = None):
    if dt is None:
        dt = state.dynamical_system.dt
    
    q_dotdots = state.get_acceleration()
    
    
    qs_next = state.qs + dt * state.q_dots
    q_dots_next = state.q_dots + dt * q_dotdots
    
    
    state_next = State(qs_next, q_dots_next)
    return state_next


def midpoint_rule_next(state, dt = None):
    if dt == None:
        dt = state.dynamical_system.dt
    
    
    
    midpoint_state = forward_euler_next(state, dt/2)
    
    
    midpoint_q_dotdots = midpoint_state.get_acceleration()
    
    qs_next = state.qs + dt * midpoint_state.q_dots
    q_dots_next = state.q_dots + dt * midpoint_q_dotdots
    
    
    state_next = State(qs_next, q_dots_next)
    return state_next



def semi_implicit_euler_next(state, dt = None):
    if dt is None:
        dt = state.dynamical_system.dt
    
    q_dotdots = state.get_acceleration()
    
    
    
    q_dots_next = state.q_dots + dt * q_dotdots
    
    qs_next = state.qs + dt * q_dots_next
    
    
    state_next = State(qs_next, q_dots_next)
    return state_next
    

def rk_slope(state):
    '''this acts as the derivative of the function y(t) = (q(t), q_dot(t)),
    which is given by f(q(t), q_dot(t)) = (q_dot(t), 1/m*forces(q,q_dot))
    '''
    f_q = state.q_dots
    f_q_dots = state.get_acceleration()
    
    return State(f_q, f_q_dots)


def rk_sampler(state, k_i, dt = None):
    '''Returns k_i_plus_one in runge kutta scheme
    
    Runge Kutta methods work by having access to the function which defines 
    the autonomous equation
        d/dt y(t) = f(y(t)).
    in this code, f is given by rk_slope.
    
    Given y(0), it then samples f at multiple points y(0) + a*dt*k_i
    to get k_i_plus_1 = f(y(0) + a*dt*k_i)
    This is the sampler that gets that value.
    
    Note how this is distinct from euler method or even midpoint method,
    and it does not help to utilize the code from forward euler.
    '''
    qs = state.qs
    q_dots = state.q_dots
    
    
    sample_qs = qs + dt * k_i.qs
    sample_q_dots = q_dots + dt * k_i.q_dots
    
    sample_state = State(sample_qs, sample_q_dots)
    
    return rk_slope(sample_state)
    

def rk_ks(state, a, dt):
    '''returns the ks for the rk-scheme with matrix a_ij
    
    The formula for the ks is
        k_l = f(y_n + dt * (a[l,0]*k_0 + a[l,1]*k_1 + dots + a[l,l-1] * k_l-1))
    an s-step scheme runs l from 0 all the way up to l=s-1. For such a scheme, a_ij
    will be an s-1 x s-1 matrix (indexed from 0 to s-2)
    '''
    s = len(a)
    k = [] #the k indexing starts at 0 instead of 1
    
    k0 = rk_slope(state)
    k.append(k0)
    #stops at s-1 as expected
    for l in range(s): 
        #stops at s-1 as expected
        
            
        sample_state_l_qs = state.qs + \
             dt * sum([a[l,m] * (k[m].qs) for m in range(l)], start=Coordinates()) 
             #sum with a start gives a 0 object to start adding to. This is especially
             #important when it comes to empty lists

        
        sample_state_l_q_dots = state.q_dots + \
            dt * sum([a[l,m] * k[m].q_dots for m in range(l)], start=Coordinates())
             
        sample_state = State(sample_state_l_qs, sample_state_l_q_dots)
        
        k_l = rk_slope(sample_state)
        k.append(k_l)
    
    return k



    

def rk_general_next(state, a, b, dt = None):
    if dt == None:
        dt = state.dynamical_system.dt
    
    s = len(b)
    
    k = rk_ks(state, a, dt)
    
    next_state_qs = state.qs + dt * sum([b[i]*k[i].qs for i in range(s)], start=Coordinates())
    next_state_q_dots = state.q_dots + dt * sum([b[i]*k[i].q_dots for i in range(s)], start=Coordinates())
    
    next_state = State(next_state_qs, next_state_q_dots)
    
    return next_state


def ssprk3_next(state, dt = None):
    if dt == None:
        dt = state.dynamical_system.dt
    
    a = np.array([[0.0, 0.0, 0.0],
                    [1.0, 0.0, 0.0],
                    [0.25, 0.25, 0.0]
                    ])
    b = np.array([1/6, 1/6, 2/3], dtype=float)

    next_state = rk_general_next(state, a, b, dt)

    return next_state


def rk45_next(state, dt = None):
    if dt == None:
        dt = state.dynamical_system.dt
    
    a = np.array([[0.0, 0.0, 0.0, 0.0],
                     [0.5, 0.0, 0.0, 0.0],
                     [0.0, 0.5, 0.0, 0.0],
                     [0.0, 0.0, 1.0, 0.0]
                     ])
    b = np.array([1/6, 1/3, 1/3, 1/6], dtype=float)

    next_state = rk_general_next(state, a, b, dt)

    return next_state



integrator_dict = {'forward euler':forward_euler_next,
               'midpoint rule': midpoint_rule_next,
               'rk45':rk45_next,
               'ssprk3':ssprk3_next,
               'verlet':verlet_next,
               'semi-implicit euler': semi_implicit_euler_next
               }



