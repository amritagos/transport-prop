# Declare exported functions
__all__ = [ 'TimeCorrelation']

import numpy as np
from ase.atoms import Atoms # ASE stuff 
import os
import io 
import re
import math

from . import misc

class TimeCorrelation():
    ''' A class that takes in two lists of Atoms objects, such that each Atom corresponds to a
    particular configuration in time. The per atom potential energies in the trajectory must be set previously. 
    By default, the number of lag times is taken to be half of the 
    steps in the trajectory, thus ensuring that every TCF calculation averages over the same number of 
    time origins. 

    traj: This must be a list of Atoms objects. Per atom energies must be present. 
    The coordinates are irrelevant from the point of view of the TCF.
    The number of frames (n_frames) is obtained from the trajectory. 

    max_tau: The maximum value of the lag time (tau). The TCF values will be averaged over different time origins, upto the max_tau.
    By default, the maximum lag time is taken to be int(0.5*n_frames). This ensures that every TCF for a particular lag time (tau) is 
    averaged over an equal number of times [equal to int(0.5*n_frames)]. TODO: allow averaging over a variable number of time origins? 

    start_t0: Indicates the first time origin from which the TCF values will be calculated. By default, this is 0, which 
    means that the first time origin will correspond to the first configuration in the trajectory

    start_tau: The first lag time to calculate; subsequent lag times are calculated to be start_tau + i * delta_tau, where i = 1, 2... 
    and delta_tau is the step size of the lag time. 

    delta_tau: The step size of the lag times. By default, this is 1. For a delta_tau of 1, the lag times will be 1, 2, 3, ..., max_tau

    '''
    def __init__(self, traj0, traj1, max_tau = None, start_t0 = 0, start_tau = 1, delta_tau = 1):
        self.traj0 = traj0 # Trajectory, list of Atoms objects ("Ground state")
        self.traj1 = traj1 # Trajectory, list of Atoms objects (so-called "excited state")
        self.n_frames = len(traj) # Number of frames in the trajectory
        self.start_t0 = start_t0 # The first time origin to calculate the MSD values from
        self.start_tau = start_tau # The first lag time 
        self.delta_tau = delta_tau # Step size in the lag times tau.  

        # TODO: error handling for start_t0 and delta_tau

        # By default, max_tau is the rounded down integer value of half of the total number of frames. 
        # Then, for every tau, the TCF is averaged over the same number of times.
        if max_tau is None:
            self.max_tau = math.floor(0.5*self.n_frames)
        elif max_tau >  math.floor(0.5*self.n_frames):
            self.max_tau = math.floor(0.5*self.n_frames) # to ensure averaging over the same number of time origins
        else:
            self.max_tau = max_tau 

        # Here, the number of time origins is equal to the maximum lag time - starting time origin index (default 0)
        self.n_origins = self.max_tau - self.start_t0
        # Assuming that the number of atoms remains constant and is the same in both trajectories
        self.n_atoms = len(self.traj0[0]) # Number of atoms 

    def e_tau_e_t0(self, f_energ_t0, f_energ_t):
        ''' This calculates the product of the fluctuations in the energy:
        f_energ_t0 (fluctuation in the energy at a particular time origin)
        and f_energ_t (at a particular lag time), 
        averaged over all the particles in the system (say, n_atoms). 
        This returns the a numPy array of size (dim, ), where dim is the number of dimensions
        i.e. 3 for a 3 D simulation.

        f_energ_t0: NumPy array of fluctuations in the energy of n_atoms at a particular time origin;
                        should have a shape of (n_atoms, dim), where dim is the number of dimensions
                        for which the TCF is being calculated.
        f_energ_t : NumPy array of fluctuations in the energy of n_atoms at time=t0+tau 
                        (where tau is a particular lag time);
                        should also have a shape of (n_atoms, dim)
        '''
        e_e0 = f_energ_t0 * f_energ_t

        # Get a single value (the numerator in the C(t) equation)
        return np.mean(e_e0, dtype=np.float64, axis=0)


    def tcf_tau(self, tau):
        ''' This calculates the TCF (aka C(t), see 10.1103/PhysRevLett.93.023004 ), for a particular lag time tau, 
        averaged over all permitted time origins [ also averaged over all the particles in the system (say, n_atoms)]. 
        At this point, we already know the first time origin from which  
        which the TCF values will be calculated from (start_t0), and we know the number of time origins (n_origins)
        This returns the averaged product of the energy fluctuations, in the same units as the per atom energy values in
        the input data.

        tau: The lag time 

        Returns a numPy array of size (1,)
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