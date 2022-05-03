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
        This returns the a numPy array of size (dim, ), where dim is the number of dimensions
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


    def vacf_tau(self, tau):
        ''' This calculates the VACF (for vx*vx, vy*vy, vz*vz), for a particular lag time tau, 
        averaged over all permitted time origins [ also averaged over all the particles in the system (say, n_atoms)]. 
        At this point, we already know the first time origin from which  
        which the VACF values will be calculated from (start_t0), and we know the number of time origins (n_origins)
        This returns the averaged squared velocity for vx*vx, vy*vy and vz*vz, in the same units as the velocities in
        the input data.

        tau: The lag time 

        Returns a numPy array of size (3,) for 3 dimensions, corresponding to vx*vx, vy*vy , vz*vz
        '''
        # VACF values of length n_origins for every time origin  
        vacf_origins = np.zeros((self.n_origins, 3)) 

        # Loop over all the time origins, starting from 
        # t0 = start_t0
        for idx in range(self.n_origins):
            t0 = idx + self.start_t0 # current time origin
            t = t0 + tau # time t = t0 + lag time 
            # Velocities at t0 and t
            vel_t0 = self.traj[t0].get_velocities() # vx, vy, vz velocities 
            vel_t = self.traj[t].get_velocities() # vx, vy, vz velocities

            # Call the VACF function for t0 and tau
            vacf_origins[idx, :] = self.v_tau_v_t0(vel_t0, vel_t)

        # Get the average over all the permitted time origins
        return np.mean(vacf_origins, dtype=np.float64, axis=0)

    def calculate_vacf(self):
        ''' Given a valid trajectory with velocities, calculate the VACF for various lag times. 
        The first lag time is given by start_tau, and subsequent values of the lag time
        are incremented in steps of delta_tau, up till the maximum lag time max_tau is reached. 

        The output numPy array will be of size (n_tau, 4), where n_tau corresponds to the number of lag times 
        calculated. The first column contains the lag times, with the other columns containing the VACF of vx*vx, vy*vy
        and vz*vz respectively.
        '''
        # where we store the lag time and the associated VACF averaged over the time origins and atoms
        vacfList = [] # will be reshaped later 
        n_tau = 0 # Number of lag times 

        # Loop over the lag times 
        for i_tau in range(self.start_tau, self.max_tau, self.delta_tau):
            n_tau += 1 
            current_tau = i_tau # Current lag time 
            # Get the MSD for this lag time over the desired time origins 
            current_vacf = self.vacf_tau(current_tau)
            # print("current tau: ", i_tau)
            # print("current vacf: ", current_vacf)
            # Update the list of VACF values (1-D for now)
            vacfList.append(current_tau)
            for k in range(3):
                vacfList.append(current_vacf[k])

        # Get a numPy array of the lag times and the VACF values
        return np.array(vacfList).reshape( (n_tau, 4) )