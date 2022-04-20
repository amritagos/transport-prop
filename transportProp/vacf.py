# Declare exported functions
__all__ = [ 'VelocityAutoCorrelation']

import numpy as np
from ase.atoms import Atoms # ASE stuff 
import os
import io 
import re
import math

from . import misc

class VelocityAutoCorrelation():
    ''' A class that takes in a list of Atoms objects, such that each Atom corresponds to a
    particular configuration in time. The velocities in the trajectory must be set previously. 
    By default, the number of lag times is taken to be half of the 
    steps in the trajectory, thus ensuring that every VACF calculation averages over the same number of 
    time origins. 

    traj: This must be a list of Atoms objects. Velocities must be present. The coordinates are irrelevant from the point of view of the VACF.
    The number of frames (n_frames) is obtained from the trajectory. 

    max_tau: The maximum value of the lag time (tau). The VACF values will be averaged over different time origins, upto the max_tau.
    By default, the maximum lag time is taken to be int(0.5*n_frames). This ensures that every VACF for a particular lag time (tau) is 
    averaged over an equal number of times [equal to int(0.5*n_frames)]. TODO: allow averaging over a variable number of time origins? 

    start_t0: Indicates the first time origin from which the VACF values will be calculated. By default, this is 0, which 
    means that the first time origin will correspond to the first configuration in the trajectory

    start_tau: The first lag time to calculate; subsequent lag times are calculated to be start_tau + i * delta_tau, where i = 1, 2... 
    and delta_tau is the step size of the lag time. 

    delta_tau: The step size of the lag times. By default, this is 1. For a delta_tau of 1, the lag times will be 1, 2, 3, ..., max_tau

    '''
    def __init__(self, traj, max_tau = None, start_t0 = 0, start_tau = 1, delta_tau = 1):
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

    def v_tau_v_t0(self, vel_t0, vel_t):
        ''' This calculates the product of velocities:
        v0 (at a particular time origin) and v_t (at a particular lag time), 
        averaged over all the particles in the system (say, n_atoms). 
        This returns the a velocity numPy array of size (dim, ), where dim is the number of dimensions
        i.e. 3 for a 3 D simulation.

        vel_t0: NumPy array of velocities of n_atoms at a particular time origin;
                        should have a shape of (n_atoms, dim), where dim is the number of dimensions
                        for which the VACF is being calculated.
        vel_t : NumPy array of velocities of n_atoms at time=t0+tau (where tau is a particular lag time);
                        should also have a shape of (n_atoms, dim)
        '''

        v_v0 = vel_t0 * vel_t

        # Get a single velocity vector vacf_vx, vacf_vy and vacf_vz
        return np.mean(v_v0, dtype=np.float64, axis=0)