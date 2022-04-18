# Declare exported functions
__all__ = [ 'MeanSquaredDisplacement']

import numpy as np
from ase.atoms import Atoms # ASE stuff 
import os
import io 
import re
import math

from . import misc

class MeanSquaredDisplacement():
    ''' A class that takes in a list of Atoms objects, such that each Atom corresponds to a
    particular configuration in time. By default, the number of lag times is taken to be half of the 
    steps in the trajectory, thus ensuring that every MSD calculation averages over the same number of 
    time origins. 

    traj: This must be a list of Atoms objects. The coordinates should be unwrapped. The number of frames (n_frames)
    is obtained from the trajectory. 

    max_tau: The maximum value of the lag time (tau). The MSD values will be averaged over different time origins, upto the max_tau.
    By default, the maximum lag time is taken to be int(0.5*n_frames). This ensures that every MSD for a particular lag time (tau) is 
    averaged over an equal number of times [equal to int(0.5*n_frames)]. TODO: allow averaging over a variable number of time origins? 

    start_t0: Indicates the first time origin from which the MSD values will be calculated. By default, this is 0, which 
    means that the first time origin will correspond to the first configuration in the trajectory

    start_tau: The first lag time to calculate; subsequent lag times are calculated to be start_tau + i * delta_tau, where i = 1, 2... 
    and delta_tau is the step size of the lag time. 

    delta_tau: The step size of the lag times. By default, this is 1. For a delta_tau of 1, the lag times will be 1, 2, 3, ..., max_tau

    dim: Indicates whether the MSD should be done for the x dimension ('x'), y dimension ('y'), z dimension ('z') or for all 
    three dimensions ('xyz'). By default, it is set to 'xyz'.
    '''
    def __init__(self, traj, max_tau = None, start_t0 = 0, start_tau = 1, delta_tau = 1, dim = 'xyz'):
        self.traj = traj # Trajectory, list of Atoms objects
        self.n_frames = len(traj) # Number of frames in the trajectory
        self.start_t0 = start_t0 # The first time origin to calculate the MSD values from
        self.start_tau = start_tau # The first lag time 
        self.delta_tau = delta_tau # Step size in the lag times tau.  

        # TODO: error handling for start_t0 and delta_tau

        # By default, max_tau is the rounded down integer value of half of the total number of frames. 
        # Then, for every tau, the MSD is averaged over the same number of times.
        if max_tau is None:
            self.max_tau = math.floor(0.5*self.n_frames)
        elif max_tau >  math.floor(0.5*self.n_frames):
            self.max_tau = math.floor(0.5*self.n_frames) # to ensure averaging over the same number of time origins
        else:
            self.max_tau = max_tau 

        # Here, the number of time origins is equal to the maximum lag time - starting time origin index (default 0)
        self.n_origins = self.max_tau - self.start_t0
        # Assuming that the number of atoms remains constant
        self.n_atoms = len(self.traj[0]) # Number of atoms 

        # To decide which dimension to calculate the MSD in: 
        if dim in ['x', 'y', 'z', 'xyz']:
            self.dim = dim
        else:
            print("You have entered an invalid option for the dimension. Defaulting to xyz.\n")
            self.dim = 'xyz'

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
        r2 = np.zeros(self.n_atoms) # Will contain the squared distance values for every atom  
    
        # Loop over all the atoms
        # to get the squared displacement given positions
        # at t0 and t
        for i in range(self.n_atoms):
            r2[i] = self.sq_disp_iatom( pos_t0[i,:], pos_t[i,:] )

        # Get the average over all atoms
        return np.mean(r2, dtype=np.float64)

    def msd_tau(self, tau):
        ''' This calculates the mean-squared displacement (MSD), for a particular lag time tau, 
        calculated for all permitted time origins, averaged over all the particles in the system (say, n_atoms). 
        At this point, we already know the first time origin from which  
        which the MSD values will be calculated from (start_t0), and we know the number of time origins (n_origins)
        This returns the squared displacement, in the same units as the coordinates in
        the input data.

        tau: The lag time 
        '''
        # MSD values of length n_origins for every time origin  
        r2_origins = np.zeros(self.n_origins) 

        # Loop over all the time origins, starting from 
        # t0 = start_t0
        for idx in range(self.n_origins):
            t0 = idx + self.start_t0 # current time origin
            t = t0 + tau # time t = t0 + lag time 
            # Positions at t0 and t
            pos_t0 = self.traj[t0].get_positions() # all x, y, z positions 
            pos_t = self.traj[t].get_positions() # all x, y, z positions
            # Get slices for dim='x', 'y', 'z': 
            pos_t0 = misc.get_slice(pos_t0, self.dim)
            pos_t = misc.get_slice(pos_t, self.dim)

            # Call the msd function for t0 and tau
            r2_origins[idx] = self.msd_tau_t0(pos_t0, pos_t)

        # Get the average over all the permitted time origins
        return np.mean(r2_origins, dtype=np.float64)

    def calculate_msd(self):
        ''' Given a valid trajectory, calculate the MSD for various lag times. 
        The first lag time is given by start_tau, and subsequent values of the lag time
        are incremented in steps of delta_tau, up till the maximum lag time max_tau is reached. 

        The output numPy array will be of size (n_tau, 2), where n_tau corresponds to the number of lag times 
        calculated. The first column contains the lag times, and the second column contains the corresponding MSD
        values. 
        '''
        # where we store the lag time and the associated MSD averaged over the time origins and atoms
        msdList = [] # will be reshaped later 

        # Loop over the lag times 
        for i_tau in range(self.start_tau, self.max_tau, self.delta_tau):
            current_tau = i_tau # Current lag time 
            # Get the MSD for this lag time over the desired time origins 
            current_msd = self.msd_tau(current_tau)
            # Update the list of MSD values (1-D for now)
            msdList.append(current_tau)
            msdList.append(current_msd)

        # Number of lag times 
        n_tau = int(0.5*len(msdList))
        # Get a numPy array of the MSD values
        msdList = np.array(msdList)
        # Reshape in the form (n_tau, 2) where every row has a
        # time lag and MSD value
        msdList = np.reshape(msdList, (n_tau,2))
        return msdList