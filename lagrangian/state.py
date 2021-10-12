#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Wed Dec 23 15:48:33 2020

@author: robertdenomme

A State object will hoold the values of the generalized coordinates, velocities 
and be responsible for computing/storing the forces 
"""



import numpy as np
import math
import copy

import itertools


class Coordinates: 
    '''use this when initializing raw data or initializing a zero coordinate state'''
    coordinates = list()
    dim = 2
    '''this should be set up by the dynamical system as a class variable'''
    
    def __init__(self, coord_dict = None):
        if coord_dict is None:
            self.coord_dict = {i:np.zeros(self.dim, dtype = float) for i in self.coordinates}
        else:
            self.coord_dict = {i:np.array(q, dtype = float) for i, q in coord_dict.items()}


    # def NewCoordinates(self, coord_dict):
    #     '''for optimization, use this method and make the qs and q_dots the correct data structure yourself'''
    #     self.coord_dict = coord_dict
        
    
    def __str__(self):
        return str(self.coord_dict)
    
    
    def copy(self):
        coord_dict_copy = copy.deepcopy(self.coord_dict)
        coordinates_copy = Coordinates(coord_dict_copy)
        
        return coordinates_copy
        
    
    
    def items(self):
        return self.coord_dict.items()
    
    
    def __iter__(self):
        return self.coord_dict.__iter__()
        
    def __getitem__(self, i):
        return self.coord_dict[i]
    
    
    def __setitem__(self, i, value):
        self.coord_dict[i] = value


    def __add__(self, other):
        coord_sum = Coordinates(self.coord_dict)
        for i in self.coordinates:
            coord_sum.coord_dict[i] += other[i]
        return coord_sum
    
    def __sub__(self, other):
        coord_sum = Coordinates(self.coord_dict)
        for i in self.coordinates:
            coord_sum.coord_dict[i] -= other[i]
        return coord_sum

    def __rmul__(self, scalar):
        coord_sum = Coordinates(self.coord_dict)
        for i in self.coordinates:
            coord_sum.coord_dict[i] *= scalar
        return coord_sum


class State:
    dynamical_system = None
    
    '''The dynamical system will be a class variable shared by all instances 
    of the class. This class variable must be set by whoever starts creating 
    these classes.
    '''
    
    def __init__(self, qs, q_dots = None):
        '''
        

        Parameters
        ----------
        qs : Dict
            coordinate dict turned into a Coordinates object which has some 
            arithmetic built in.
        q_dots : same dict type as qs, optional
            DESCRIPTION. The default is None.

        Returns
        -------
        None.

        '''
        self.qs = Coordinates(qs)
        self.q_dots = Coordinates(q_dots)
        
        forces = dict()#? property
    
    
    # def NewState(self, qs, q_dots):
    #     '''for optimization, use this method and make the qs and q_dots the correct data structure yourself'''
    #     self.qs = Coordinates.NewCoordinates(qs)
    #     self.q_dots = Coordinates.NewCoordinates(q_dots)
    #     return self
    
    
    @property
    def masses(self):
        return self.dynamical_system.masses
    
    
    def copy(self):
        qs_copy = self.qs.copy()
        q_dots_copy = self.q_dots.copy()
        state_copy = State(qs_copy, q_dots_copy)
        return state_copy
        
    
    
    def get_forces(self):
        state = self
        dynamical_system = self.dynamical_system 
        
        
        forces = Coordinates()
        
        
        
        '''Compute the force generated by the potentials'''
        for _, potential in dynamical_system.potentials.items():
            for i in potential.coordinates:
                force = - potential.gradient[i](state.qs)
                forces[i] += force
        
        if dynamical_system.cell_list_potentials:
            '''Compute the forces generated by cell-list potentials next'''
            forward_neighbor_lists = dynamical_system.forward_neighbor_lists
            for q in state.qs:
                for p in forward_neighbor_lists[q]:
                    for potential in dynamical_system.cell_list_potentials.values():
                        force = - potential.gradient[0](state.qs[q], state.qs[p])
                        #need to define a cell_list_pair_potential object...
                        forces[q] += force
                        forces[p] -= force
        
        '''Compute the force generated by the constraints'''
        
        
        masses = self.masses
        #masses indexed by coordinates in dynamical_system.coordinates
        
        
        #need an enumeration for the constraints, indexed by k for the equation
        #and j for the lagrange multiplier
        N_constraints = len(dynamical_system.constraints)
        
        
        '''Start with the RHS. We'll leave it as a dictionary type 
        object at first and then convert to a matrix later
        '''
        RHS = np.zeros(N_constraints)
        for k, F_k in enumerate(dynamical_system.constraints):
            '''term 1.'''
            hessian_k = F_k.get_hessian(state.qs)
            RHS[k] =  - hessian_k.double_dot(state.q_dots)
            
            '''term 2.'''
            for potential in dynamical_system.potentials:
                for i in potential.coordinates():
                    
                    RHS[k]  = RHS[k] - \
                        masses[i]**-1 * \
                            np.dot( \
                                   F_k.gradient[i](state.qs), \
                                       potential.gradient[i](state.qs)
                                       )
                            
            
        '''Now the LHS. This will be a double dict type object 
        for rows indexed by constraint equations and columns indexed by 
        lagrange multipliers'''
        
        Matrix = np.zeros((N_constraints,N_constraints))
        for k, F_k in enumerate(dynamical_system.constraints):
            for j, F_j in enumerate(dynamical_system.constraints):
                
                for i in F_k.coordinates():
                    #note the row-column indexing k,j
                    Matrix[k,j] += masses[i]**-1 * \
                        np.dot(\
                               F_k.gradient[i](state.qs),
                               F_j.gradient[i](state.qs)
                               )
            
        Lambda = np.linalg.solve(Matrix, RHS)
        
        for k, F_k in enumerate(dynamical_system.constraints):
            for i in dynamical_system.coordinates:
                forces[i] += Lambda[k] * F_k.gradient[i](state.qs)
                
        return forces
            
    def get_acceleration(self):
        forces = self.get_forces()
        masses = self.masses
        accel_qs = Coordinates(
            {i:masses[i]**-1 * forces[i] for i in forces.coordinates}
            )
        return accel_qs


class TrajectoryData:
    '''Object that holds the states that occur during a simulation over time. Also used for rendering
    
    For each time step a trajectory should hold a state of the system.
    It will also hold list of particles, 
    their types an optional color for each type and the coordinates for each
    particle. 
    
    It should have methods for getting the various data used in rendering.
    needed 
    '''
    def __init__(self, dynamical_system):
        self.dynamical_system = dynamical_system
        self.states = [dynamical_system.initial_state]
        
        self.dt = dynamical_system.dt
        
        
        self.FRAME_RATE = 1/20
        self.skipframes = math.ceil(self.FRAME_RATE / self.dt)
        self.frame_dt = self.dt * self.skipframes
    
    
    def append(self, value):
        self.states.append(value)
    
    
    def __getitem__(self, i):
        return self.states[i]
    
    def process_data_for_rendering(self):
        '''Run this before heading to the renderer to assemble data ahead of time'''
        
        self.N_time_steps = len(self.states)
        self.N_frames = len(list(range(0, self.N_time_steps, self.skipframes)))
        N_frames = self.N_frames
        
        
        
        coordinates = self.dynamical_system.coordinates
        
        
        
        vertices_at_time = []
        for t, time_step in zip(range(0,self.N_time_steps,self.skipframes), itertools.count(start=0)):
            vertices_at_time.append([]) # adding t'th component to thislist of lists
            for i in self.dynamical_system.rendered_paths:
                vertices_at_time[time_step].append(\
                                           (self.states[t].qs[i][0],
                                            self.states[t].qs[i][1])
                                           )
                
        self.vertices_at_time = vertices_at_time
        
        
        # self.xs_at_time = [[self.states[t].qs[i][0] for i in coordinates] for t in range(0, self.N_time_steps, self.skipframes)]
        # self.ys_at_time = [[self.states[t].qs[i][1] for i in coordinates] for t in range(0, self.N_time_steps, self.skipframes)]
        
        
    def get_xs_at_time(self, time):
        '''Get all x-values at given time'''
        return self.xs_at_time[time]
    
    
    def get_ys_at_time(self, time):
        '''Get all x-values at given time'''
        return self.ys_at_time[time]
   
    def get_x(self, i, t):
        return self.states[t].qs[i][0]
    
    def get_y(self, i, t):
        return self.states[t].qs[i][1]