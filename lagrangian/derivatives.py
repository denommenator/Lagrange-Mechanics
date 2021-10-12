#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Sun Dec  6 12:59:55 2020

@author: robertdenomme

Derivatives module. See test file 
moleculardynamics.tests.numericaldifferentiation for more details

basic usage:
    
f = lambda x,y,z, G=4: x*y*z*G**2
f = derivatives.add_partials(f)

f.partial[0](1,1,1,G=4)
f.partial['G'](1,1,1,G=4)
"""

import numpy as np
import math

import functools


import inspect

np.seterr(divide = 'raise', invalid = 'raise')


def single_variable_difference_quotient(f, x):
        """See wikipedia article on numerical differentiation for this choice 
        of dx. Basically we use square root of machine epsilon
        since you can safeuly multiply this by anything, then
        get a dx based on it scaled by your x since floating point bs.
        
        
        """
        
        epsilon = np.finfo(float).eps
        h = math.sqrt(epsilon) * x
        
        xph = x + h
        xmh = x - h
        
        two_h = xph - xmh
        try:
            return (f(xph)-f(xmh))/two_h
        
        except (FloatingPointError, ZeroDivisionError):
            h = math.sqrt(epsilon)
            return (f(h)-f(epsilon)) / (h-epsilon)

def partial_keyword(key, f):
    """returns partial derivative of f w.r.t. a keyword argument
    modeled on 
    single_variable_difference_quotient above"""
    def F_key (*args, **kwargs):
        arg_key = kwargs[key]
        
        
        def f_args_held_constant(x): 
            kwargs[key]=x
            return f(*args, **kwargs)
        #note the use of 'closures' in the previous line. Outside the function 
        #partial i,f the other arguments will still be well defined.
        
        return single_variable_difference_quotient(f_args_held_constant, arg_key)
    return F_key


def partial(i,f):
    """returns partial derivative of f modeled on 
    single_variable_difference_quotient above"""
    def F_i (*args, **kwargs): 
        pre_args = args[:i]
        arg_i = args[i]
        post_args = args[i+1:]
        
        
        f_args_held_constant = lambda x : f(*pre_args, x , *post_args, **kwargs)
        #note the use of 'closures' in the previous line. Outside the function 
        #partial i,f the other arguments will still be well defined.
        
        return single_variable_difference_quotient(f_args_held_constant, arg_i)
    return F_i
        

class Partials:
    def __init__(self,f):
        self.f = f
        
        
    def __getitem__(self, index):
        if  isinstance(index, int):
            return partial(index, self.f)
            
        else:
            return partial_keyword(index, self.f)
            
        
def add_partials(f):
    """This function returns amodified version of the function f, 
    where calling f(x) results in the normal result, but f.partial[i](x)
    returns the partial derivative of f w.r.t. the i-th positional variable
    of f evaluated at x. 
    Note that add_partials is idempotent, and that if a user has already 
    defined explicit partial derivatives (which might evaluate slightly faster?)
    then it will leave them alone."""
    try: 
        """check if f already comes with its partial derivatives defined."""
        f.partial
        
    
    
    except AttributeError:
        """partials are not defined, so let's calculate/store them
        using difference quotients"""
        
        f.partial = Partials(f)
    finally:
         return f
    

'''The next section defines higher order gradients of functions w.r.t.
variables q which themselves have multiple variables q = (x,y,z)
'''



# class Gradients:
#     def __init__(self,f):
#         self.f = f
        
        
#     def __getitem__(self, index):
#         if  isinstance(index, int):
#             return partial(index, self.f)
            
#         else:
#             return partial_keyword(index, self.f)
       



def two_d_gradient(f, q):
    '''Returns the value of gradient of the function f at the point q.
    
    q is assumed to be a 2d vector
    '''
    q_x = q[0]
    q_y = q[1]
    
    f_restricted_x = lambda x: f(np.array([x, q_y]))
    f_restricted_y = lambda y: f(np.array([q_x, y]))
    
    f_x = single_variable_difference_quotient(f_restricted_x, q_x)
    f_y = single_variable_difference_quotient(f_restricted_y, q_y)
    
    
    return np.array([f_x, f_y])
    

def get_gradient_i(U, i):
    def U_gradient_i(*args, **kwargs):
        '''returns a function which gives the gradient of U wrt. qi at the input determined by *args, **kwargs
        '''
        def U_restricted_i(*args, **kwargs):
            '''Returns a function of a signel variable q=qi (ignoring the value of qi passed in)
            
            Restricts U to a function of the variable, qi while holding the other
            values of qj passed in constant
            '''
            pre_args = args[:i]
            arg_i = args[i] #actually a throw-away variable
            post_args = args[i+1:]
            
            U_args_held_constant = lambda qi : U(*pre_args, qi , *post_args, **kwargs)
            
            return U_args_held_constant
        
        
        
        
        U_restricted = U_restricted_i(*args, **kwargs)
        
        return two_d_gradient(U_restricted, args[i]) 
        
    return U_gradient_i
        
        

def get_gradient_functions(U):
    ''' Assuming U is a function of coordinates q1, q2, this function adds gradients wrt these coords
    
    The input to U should be either a list, or a numpy array of coordinates q1, ... qn
    each with dimension 2. The gradient U.gradient[i] should return a 2d numpy array
    partial qi_x, partial qi_y.
    
    '''
    
    sig = inspect.signature(U)
    n =len([param for param in sig.parameters.values() 
            if (param.kind == param.POSITIONAL_ONLY or \
                param.kind == param.POSITIONAL_OR_KEYWORD)
                ])
    '''n is the number of non-keyword only arguments, and we will define
    the gradient wrt any non-keyword only argument. they are only executed 
    if asked for via U.gradient[i](q1,q2,...) so there is no trouble
    if you have some positional or keyword arguments you don't want the gradient wrt.
    '''
    gradients = []
    
    for i in range(n):
        gradients.append(get_gradient_i(U,i))
            
    
    return gradients



def add_gradients(U):
    """ after calling U = add_gradients(U), you'll be able to use U.gradient[i](q1,q2,...)
    
    This function returns a modified version of the function U, 
    where calling U(q1,q2,...) results in the normal result, but U.gradient[i](q1,q2,...)
    returns the gradient of U w.r.t. the i-th positional variable
    of U evaluated at (q1,q2,...). 
    Note that add_gradients is idempotent, and that if a user has already 
    defined explicit gradients (which might evaluate slightly faster?)
    then it will leave them alone."""
    try: 
        """check if f already comes with its partial derivatives defined."""
        U.gradient
        
    
    
    except AttributeError:
        """partials are not defined, so let's calculate/store them
        using difference quotients"""
        
        U.gradient = get_gradient_functions(U)
    finally:
         return U


    
    
    #you can actually get the name of the function using this code
    #return [f'{f.__name__}_0', f'{f.__name__}_1', f'{f.__name__}_1', '...']