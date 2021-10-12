#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Mon Dec 21 19:38:13 2020

@author: robertdenomme

A Potential object represents a potential energy that depends on some 
collection of coordinates. It has an attribute .gradient which you can use
as potential.gradient[i](state)
"""
import numpy as np
import math

import lagrangian.derivatives as derivatives

np.seterr(divide = 'raise', invalid = 'raise')

epsilon = np.finfo(float).eps

def single_variable_difference_quotient(f, x):
        """See wikipedia article on numerical differentiation for this choice 
        of dx. Basically we use square root of machine epsilon
        since you can safeuly multiply this by anything, then
        get a dx based on it scaled by your x since floating point bs.
        
        
        """
        
        h = math.sqrt(epsilon) * x
        
        xph = x + h
        xmh = x - h
        
        two_h = xph - xmh
        try:
            return (f(xph)-f(xmh))/two_h
        
        except (FloatingPointError, ZeroDivisionError):
            h = math.sqrt(epsilon)
            return (f(h)-f(epsilon)) / (h-epsilon)


def single_var_gradient(f, q):
    '''Returns the value of gradient of the function f at the point q.
    
    q is assumed to be a 2d vector
    '''
    q_x = q[0]
    q_y = q[1]
    
    def f_restricted_x(x): return f(np.array([x, q_y]))
    def f_restricted_y(y): return f(np.array([q_x, y]))
    
    f_x = single_variable_difference_quotient(f_restricted_x, q_x)
    f_y = single_variable_difference_quotient(f_restricted_y, q_y)
    
    
    return np.array([f_x, f_y])
    
    

class Gradient:
    
    
    def __init__(self, function, args_list, **kwargs):
        self.function = function
        self.args_list = args_list
        self.kwargs = kwargs
    
    
    def __getitem__(self, i):
        
        def grad_i(qs):
            qs = qs
            args = {i : qs[i] for i in self.args_list}
            kwargs = self.kwargs
            function = self.function
            
            def single_variable_function(x) :
                args[i] = x
                args_tuple = (args[i] for i in args.keys())
                return function(*args_tuple, **kwargs)
            
            return single_var_gradient(single_variable_function, qs[i])
            
        
        return grad_i 
    
    
    def add_gradients(self, function, args_list, **kwargs):
        function.gradient = self.Gradient(function, args_list, **kwargs)
        return function
    


class Potential:
    def __init__(self, potential_function, coordinates, **kwargs):
        '''The args are strings that hold the name of each coordinate appearing in the input variables to the potential function.
        
        
        **kwargs should hold any additional arguments passed to the potential.
        We might need to manage where they come from...'''
        
        self.potential_function = potential_function
        self.coordinates = coordinates
        self.kwargs = kwargs
        
        self.gradient = Gradient(potential_function, coordinates, **kwargs)
        
        
        
        
        
    def __call__(self, qs):
        args = (qs[i] for i in self.coordinates)
        return self.potential_function(*args, **self.kwargs)
    
    
    
class CellListsPairPotential:
    def __init__(self, potential_function, **kwargs):
        '''potential function is just a normal function of two variables.
        
        Should be symmetric between the two variables.
        
        **kwargs should hold any additional arguments passed to the potential.
        We might need to manage where they come from...'''
        
        self.potential_function = derivatives.add_gradients(potential_function)
        self.kwargs = kwargs
        self.gradient = self.potential_function.gradient
        
        
     
    pass
    
class Constraint:
    def __init__(self, constraint_function, coordinates, **kwargs):
        '''The args are strings that hold the name of each coordinate appearing in the input variables to the potential function.
        
        
        **kwargs should hold any additional arguments passed to the potential.
        We might need to manage where they come from...'''
        
        self.constraint_function = constraint_function
        self.coordinates = coordinates
        self.kwargs = kwargs
        
        self.gradient = Gradient(constraint_function, coordinates, **kwargs)
        
        
        
        
        
    def __call__(self, qs):
        args = (qs[i] for i in self.coordinates)
        return self.constraint_function(*args, **self.kwargs)
    

    