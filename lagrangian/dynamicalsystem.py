#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Mon Dec 21 17:50:33 2020

@author: robertdenomme


The dynamical system holds the information about the system that does
not depend on time. 
"""
import math
from tqdm import tqdm

from matplotlib.path import Path


import lagrangian.potentials

from lagrangian.state import State, Coordinates, TrajectoryData
import lagrangian.collisions as collisions
import lagrangian.integrators as integrators
from lagrangian.renderer import MatPlotRenderer
import lagrangian.celllists as celllists

integrator_dict = integrators.integrator_dict 

class DynamicalSystem:
    
    
    def __init__(self, 
                 qs_init_raw, 
                 q_dots_init_raw = None, 
                 masses = None,
                 dt = 1/60,
                 xlim = (-10, 10),
                 ylim = (-10, 10),
                 cell_list_dx = 1,
                 cell_list_dy = 1,
                 wall_elasticity = 1,
                 integrator_code = 'ssprk3'):
        
        
        
        self.qs_init_raw = qs_init_raw
        self.q_dots_init_raw = q_dots_init_raw
        self.masses = masses or {i:1 for i in qs_init_raw.keys()}
        
        self.xlim = xlim
        self.ylim = ylim
        self.wall_elasticity = wall_elasticity
        
        self.cell_list_dx = cell_list_dx
        self.cell_list_dy = cell_list_dy
        
        
        
        '''setup class variables in state and coordinate modules'''
        self.coordinates = list(self.qs_init_raw.keys())
        
        self.rendered_paths = []
        self.rendered_path_codes = []
        '''rendered_paths is a list of the paths you would like to draw. 
        It has two dimensions in case you want to draw disconnected paths:
        
        -paths - used to hold the various paths you want rendered
        -path node list
        
        Example 1: To render 2 line segments between 2 sets of particles, q1, q2, p1, p2, p3:
            rendered_paths =[#two paths
                             ['q1','q2'],
                             ['p1','p2', 'p3']
                             ]
                                  
        Example 2: To render a cycle of coordinates q1,q2,q3,...
            rendered_paths = [#one path
                               ['q1','q2','q3',..., 'qN','q1']
                               ]
            '''
        
        
        
        Coordinates.coordinates = self.coordinates
        State.dynamical_system = self
        
        
        self.initial_state = State(qs_init_raw, q_dots_init_raw)
        
        ''''set up instance variables for this dynamical system'''
        self.potentials = dict()
        self.cell_list_potentials = dict()
        self.constraints = dict()
        self.dt = dt
        
        try:
            self.integrator = integrator_dict[integrator_code]
        except(KeyError):
                raise KeyError(integrator_code, 
                               ' not found. The valid integrator keys are ',
                               integrator_dict.keys())
        
        
        
        
        
    def add_potential(self, potential_function, args_list, **kwargs):
        potential = lagrangian.potentials.Potential(potential_function, args_list, **kwargs)
        potential_index = len(self.potentials)
        self.potentials[potential_index] = potential
    
    
    def add_cell_list_pair_potential(self, potential_function, **kwargs):
        potential = lagrangian.potentials.CellListsPairPotential(potential_function, **kwargs)
        potential_index = len(self.cell_list_potentials)
        self.cell_list_potentials[potential_index] = potential
        
    
    def add_constraint(self, constraint_function, args_list, **kwargs):
        constraint = lagrangian.potentials.Constraint(constraint_function, args_list, **kwargs)
        constraint_index = len(self.constraints)
        self.constraints[constraint_index] = constraint
    
    
    
    
    def test(self):
        self.initial_state = State(self.qs_init_raw, self.q_dots_init_raw)
        self.initial_state.compute_forces()
        return self.initial_state.forces
    
    
    
    def add_rendered_path(self, path_list):
        '''must pass the path list as an actual list'''
        self.rendered_paths += path_list
        self.rendered_paths.append(path_list[0]) #the extra vertex is a dummy one used for the stop code
        
        self.rendered_path_codes.append(Path.MOVETO)
        self.rendered_path_codes += [Path.LINETO] * (len(path_list) - 1)
        self.rendered_path_codes.append(Path.CLOSEPOLY)
    
        
    def run_dynamics(self, total_time):
        '''Iterate the system forward in time N_steps time steps and store in a TrajectoryData object
        '''
        if self.cell_list_potentials:
            self.run_cell_list_dynamics(total_time)
        else:
            
            dt = self.dt
            N_steps = math.ceil(total_time / dt )
            
        
            self.trajectory_data = TrajectoryData(self)
            t=0 + dt
            second_state = integrators.midpoint_rule_next(self.trajectory_data[0], dt)
            self.trajectory_data.append(second_state)
            t += dt
            
            for i in tqdm(range(1, N_steps-1)):  #add tqdm loop counter display
                next_state = self.integrator(self.trajectory_data[-1], dt)
                
                next_state = collisions.resolve_wall_collisions(next_state,
                                                           xlim = self.xlim,
                                                           ylim = self.ylim,
                                                           wall_elasticity = self.wall_elasticity)
                self.trajectory_data.append(next_state)
                
                t += dt
                
        print("dynamics finished!\r")
            
        #self.trajectory_data = trajectories
    
    def run_cell_list_dynamics(self, total_time):
        dt = self.dt
        N_steps = math.ceil(total_time / dt )
        
    
        self.trajectory_data = TrajectoryData(self)
        t=0 + dt
        
        cell_lists = celllists.CellLists(self.xlim, 
                                         self.cell_list_dx,
                                         self.ylim,
                                         self.cell_list_dy)
        self.forward_neighbor_lists = \
            cell_lists.generate_forward_neighbor_lists(self.trajectory_data[-1])
        
        
        second_state = integrators.midpoint_rule_next(self.trajectory_data[0], dt)
        self.trajectory_data.append(second_state)
        t += dt
        
        for i in tqdm(range(2, N_steps-1)):  #add tqdm loop counter display
            if i % 5 == 0:
                self.forward_neighbor_lists = \
                    cell_lists.generate_forward_neighbor_lists(self.trajectory_data[-1])
            next_state = self.integrator(self.trajectory_data[-1], dt)
            #print(next_state.qs)
            
            next_state = collisions.resolve_wall_collisions(next_state,
                                                       xlim = self.xlim,
                                                       ylim = self.ylim,
                                                       wall_elasticity = self.wall_elasticity)
            self.trajectory_data.append(next_state)
            
            t += dt
    
    def display(self):
        trajectories = self.trajectory_data
        
        renderer = MatPlotRenderer(trajectories)

        renderer.display()
    