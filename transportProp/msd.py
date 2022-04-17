# Declare exported functions
__all__ = [ 'MeanSquaredDisplacement']

import numpy as np
from ase.atoms import Atoms # ASE stuff 
import os
import io 
import re
import math

class MeanSquaredDisplacement():
    ''' A class that takes in a list of Atoms objects, such that each Atom corresponds to a
    particular configuration in time. By default, the number of lag times is taken to be half of the 
    steps in the trajectory, thus ensuring that every MSD calculation averages over the same number of 
    time origins. 

    traj: This must be a list of Atoms objects. The coordinates should be unwrapped. The number of frames (n_frames)
    is obtained from the trajectory. 

    max_tau: The maximum value of the lag time (tau). The MSD values will be averaged over different time origins, upto the max_tau.
    By default, the maximum lag time is taken to be int(0.5*n_frames). This ensures that every MSD for a particular lag time (tau) is 
    averaged over an equal number of times [equal to int(0.5*n_frames)]

    start_tau: Indicates the first time origin from which the MSD values will be calculated. By default, this is 0, which 
    means that the first time origin will correspond to the first configuration in the trajectory

    delta_tau: The step size of the lag times. By default, this is 1. For a delta_tau of 1, the lag times will be 1, 2, 3, ..., max_tau
    '''
    def __init__(self, traj, max_tau = None, start_tau = 0, delta_tau = 1):
        self.traj = traj # Trajectory, list of Atoms objects
        self.n_frames = len(traj) # Number of frames in the trajectory
        self.start_tau = start_tau # The first time origin to calculate the MSD values from
        self.delta_tau = delta_tau # Step size in the lag times tau.  

        # By default, max_tau is the rounded down integer value of half of the total number of frames. 
        # Then, for every tau, the MSD is averaged over the same number of times.
        if max_tau is None:
            max_tau = math.floor(0.5*self.n_frames)          

    def sq_disp_iatom(self, pos_t0_iatom, pos_t_iatom):
        ''' Returns the squared displacement (in the same units as the coordinates) 
        given the (unwrapped) positions of a particular atom, at time=t0 (the time origin)
        and time=t0+tau (where tau is a particular lag time). 

        pos_t0_iatom: Unwrapped coordinates for a particular atom, which  
                            should be of size (dim, ), where dim is the number of dimensions. 
        pos_t_iatom: Unwrapped coordinates for the atom iatom, which should be
                            of size (dim, ).
        '''
        dr = pos_t0_iatom - pos_t_iatom # (dim, ) difference in unwrapped coord for iatom
        r = np.linalg.norm(dr)
        return r**2 # return MSD_i between t0 and t

    def msd_tau_t0(self, pos_t0, pos_t):
        ''' This calculates the mean-squared displacement (MSD), given a particular time origin,
        and a particular lag time, averaged over all the particles in the system (say, n_atoms). 
        This returns the squared displacement, in the same units as the coordinates in
        the input data.

        pos_t0: NumPy array of positions of n_atoms at a particular time origin;
                        should have a shape of (n_atoms, dim), where dim is the number of dimensions
                        for which the MSD is being calculated.
        pos_t : NumPy array of positions of n_atoms at time=t0+tau (where tau is a particular lag time);
                        should also have a shape of (n_atoms, dim)
        '''
        n_atoms = np.shape(pos_t0)[0] # Number of total atoms 
        r2 = np.zeros(n_atoms) # Will contain the squared distance values for every atom  
    
        # Loop over all the atoms
        # to get the squared displacement given positions
        # at t0 and t
        for i in range(n_atoms):
            r2[i] = self.sq_disp_iatom( pos_t0[i,:], pos_t[i,:] )

        # Get the average over all atoms
        return np.mean(r2, dtype=np.float64)