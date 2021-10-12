#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Thu Dec 31 12:52:53 2020

@author: robertdenomme

Cell Lists version 2.0 attempts to actually make an O(N) implementation
of cell lists

"""



import math
from itertools import product
import numpy as np


from lagrangian.state import State


def forward_neighbor_cells(i, j):
        return set([(i+1, j-1), 
                   (i+1, j),
                   (i+1, j+1),
                   (i, j+1)])


class CellLists:
    '''CellLists partitions region into cells and computes neighbor cells
    
    CellLists partitions region xlim x ylim into cells indexed by (i,j).
    Specifically, cell (i,j) is given by,
        [xlim[0]+i*dx, xlim[0]+(i+1)*dx) x [ylim[0]+j*dy, ylim[0]+(j+1)*dy)
    The two lists give coordinates for the qs in each cell
        
        qs_cell[(i,j)]
        qs_forward_neighbor_cells[(i,j)] 
        
    '''
    
    def __init__(self, xlim, suggested_dx, ylim, suggested_dy):
        '''Initialize and start computing the cell lists. 
        
        Note that dx and 
        dy might be raised slightly in order to make the cells divide
        the limits exactly.'''
        self.xlim = xlim
        self.ylim = ylim
        
        
        
        
        self.dx, self.N_x, self.dy, self.N_y = \
            self.determine_cell_dims(suggested_dx, suggested_dy)
        '''Computes and defines the following attributes:
                dx, N_x, dy, N_y
        '''
        
        self.cell_lists = dict()
        
        
        
        
    def determine_cell_dims(self, suggested_dx, suggested_dy):
        N_x_inside_boundary = math.floor((self.xlim[1]-self.xlim[0]) / suggested_dx) 
        #boundary refers
        dx = (self.xlim[1]-self.xlim[0]) / N_x_inside_boundary
        
        N_y_inside_boundary = math.floor((self.ylim[1]-self.ylim[0]) / suggested_dy)
        dy = (self.ylim[1]-self.ylim[0]) / N_y_inside_boundary 
        
        N_x, N_y = N_x_inside_boundary + 2, N_y_inside_boundary + 2
        #because we are going to add a single unbouded cell containing all x's to the left of xlim[0]
        # and similarly for above x[1] and for the y's.
        return dx, N_x, dy, N_y
    
    
    
        
    def cell_of_q(self, q):
        '''returns the index (i,j) of the cell q = np.array of dim 2 resides in'''
        xlim, dx, ylim, dy = self.xlim, self.dx, self.ylim, self.dy
        cell_x = math.floor((q[0]-xlim[0])/dx) + 1
        ''''the following logic speels out what happens if the q is outside the defined region'''
        if cell_x <= 0:
            cell_x = 0
        elif cell_x >= self.N_x-1:
            cell_x = self.N_x-1
            
        cell_y = math.floor((q[1]-ylim[0])/dy) + 1
        if cell_y <= 0:
            cell_y = 0
        elif cell_y >= self.N_y-1:
            cell_y = self.N_y-1
        return (cell_x, cell_y)
    
    
    def _generate_cell_lists(self, state):
        '''This function generates the cell lists cell_lists[(i,j)] 
         
        cell_lists[(i,j)] is a a dictionary which holds the names of the 
        qs in cell i,j.
        Cells with no qs in them will have no key (i,j) in this dictionary.
        
        This function returns the cell_lists object, but also stores it for 
        future use in the cell lists object
        '''
        self.cell_lists = dict() #make sure to clear the old cell lists
        for q in state.qs.coordinates:
            (i,j) = self.cell_of_q(state.qs[q])
            if (i,j) in self.cell_lists:
                self.cell_lists[(i,j)].append(q)
            else:
                self.cell_lists[(i,j)] = [q]
        
        return self.cell_lists
    
    
    
    
    
    def generate_forward_neighbor_lists(self, state):
        '''This function generates the forward_neighbor_list[q]
        
        forward_neighbor_list is a dictionary, forward_neighbor_list[q]
        holds all the forward neighbors of q.
        
        This function returns the forward_neighbor_lists object, but also stores it for 
        future use in the cell lists object
        '''
        self._generate_cell_lists(state)
        self.forward_neighbor_lists = dict() #make sure to clear the old cell lists
        
        for q in state.qs.coordinates:
            self.forward_neighbor_lists[q] = []
            (i,j) = self.cell_of_q(state.qs[q])
            for (k,l) in forward_neighbor_cells(i,j):
                if (k,l) in self.cell_lists:
                    self.forward_neighbor_lists[q] += self.cell_lists[(k,l)]
            for p in self.cell_lists[(i,j)]:
                if p>q:
                    self.forward_neighbor_lists[q].append(p)
        return self.forward_neighbor_lists
        
        
if __name__=='__main__':
    
    qs = [[0.5,4.1], 
          [2.3, 5.6], 
          [2.4, 6.5],
          [2.5, 4.5],
          [3.5, 4.5],
          [3.5, 5.5],
          [3.5, 6.5],
          [4.5, 3.6], 
          [5.9, 5.4], 
          [6.01, 1.5]]
    
    qs = {'q_'+str(i):q_i for i, q_i in enumerate(qs)}
    
    
    
    initial_state = State(qs)
    initial_state.qs.coordinates = initial_state.qs.coord_dict.keys()
    
    cell_lists = CellLists(
                           xlim=(1,6), 
                           suggested_dx=1, 
                           ylim=(1,6), 
                           suggested_dy=1)
    
    cell_lists.generate_cell_lists(initial_state)
    cell_lists.generate_forward_neighbor_lists(initial_state)
    
    print(cell_lists.forward_neighbor_lists)
    
    for i in qs.keys():
        print(i, ': ', cell_lists.cell_of_q(qs[i]), 'neighbors: ', 
              cell_lists.forward_neighbor_lists[i])
    
    
    
        
            
        