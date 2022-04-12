# Declare exported functions
__all__ = [ 'tau_t0', 'my_read_lammps_data']

import numpy as np
from ase.atoms import Atoms # ASE stuff 
import os
import io 
import re
import math

def sq_disp_iatom(pos_t0_iatom, pos_t_iatom):
    ''' Returns the squared displacement (in the same units as the coordinates) 
    given the (unwrapped) positions of a particular atom, at time=t0 (the time origin)
    and time=t0+tau (where tau is a particular lag time). 

    pos_t0_iatom: Unwrapped coordinates for a particular atom, which  
                            should be of size (dim, ), where dim is the number of dimensions. 
    pos_t_iatom: Unwrapped coordinates for the atom iatom, which should be
                            of size (dim, ).
    '''

def tau_t0(pos_t0, pos_t):
    ''' This calculates the mean-squared displacement (MSD), given a particular time origin,
    and a particular lag time, averaged over all the particles in the system (say, n_atoms). 
    This returns the squared displacement, in the same units as the coordinates in
    the input ASE Atoms objects.

    pos_t0: NumPy array of positions of n_atoms at a particular time origin;
                    should have a shape of (n_atoms, dim), where dim is the number of dimensions
                    for which the MSD is being calculated.
    pos_t : ASE Atoms object at time=t0+tau (where tau is a particular time lag);
                    should also have a shape of (n_atoms, dim)
    '''
    n_atoms = np.shape(pos_t0)[0] # Number of total atoms 