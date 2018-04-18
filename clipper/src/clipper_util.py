# @Author: Tristan Croll
# @Date:   16-Oct-2017
# @Email:  tic20@cam.ac.uk
# @Last modified by:   Tristan Croll
# @Last modified time: 18-Apr-2018
# @License: Creative Commons BY-NC-SA 3.0, https://creativecommons.org/licenses/by-nc-sa/3.0/.
# @Copyright: Copyright 2017-2018 Tristan Croll



import numpy
from clipper_python import Util as _util

### Constants ###
eightpi2 = _util.eightpi2()

### Math functions ###

def atanh(x):
    ''' Hyperbolic tan '''
    return _util.atanh(x)

def b2u(b_factor):
    ''' Convert isotropic B-factor to U-value '''
    return _util.b2u(b_factor)

def u2b(u_value):
    ''' Convert isotropic U-value to B-factor '''
    return _util.u2b(u_value)

def bessel_i0(x):
    ''' Modified Bessel function of the first kind '''
    return _util.bessel_i0(x)

def d2rad(angle):
    ''' Convert degrees to radians '''
    return _util.d2rad(angle)

def rad2d(angle):
    ''' Convert radians to degrees '''
    return _util.rad2d(angle)

def int_ceiling(fp_num):
    ''' Round the given number up and return as an integer. '''
    return _util.intc(fp_num)

def int_floor(fp_num):
    ''' Round the given number down and return as an integer. '''
    return _util.intf(fp_num)

def int_round(fp_num):
    ''' Round the given number to the nearest integer. '''
    return _util.intr(fp_num)

def invsim(x):
    ''' Inverse Sim function: I1(X)/I0(X). '''
    return _util.invsim(x)



### Useful utility functions ###
def get_minmax_grid(coords, cell, grid_sampling):
    '''
    Return a numpy array containing the (min, max) grid coordinates
    of a box encompassing the given coordinates.
    Args:
        coords ( [n*3 float]):
            An n*3 array of coordinates
        cell (clipper.Cell object)
        grid_sampling (clipper.Grid_sampling object)
    '''
    return _util.get_minmax_grid(coords, cell, grid_sampling)
