# Declare exported functions
__all__ = [ 'SolvationTimeCorrelation']

import numpy as np
from ase.atoms import Atoms # ASE stuff 
import os
import io 
import re
import math
import numpy as np

from . import misc

class SolvationTimeCorrelation():
    ''' A class that takes in a numPy array of potential energy per atom differences (n_atoms, n_frames), 
    such that each column corresponds to a particular configuration in time. 
    By default, the number of lag times is taken to be half of the 
    steps in the trajectory, thus ensuring that every TCF calculation averages over the same number of 
    time origins. 

    energy_diff_array: This is a numPy array with the energy differences for each atom (n_atoms, n_frames)
    The coordinates are irrelevant from the point of view of the TCF.
    The number of frames (n_frames) is obtained from the number of columns (counting from start_t0). 

    max_tau: The maximum value of the lag time (tau). The TCF values will be averaged over different time origins, upto the max_tau.
    By default, the maximum lag time is taken to be int(0.5*n_frames). This ensures that every TCF for a particular lag time (tau) is 
    averaged over an equal number of times [equal to int(0.5*n_frames)]. TODO: allow averaging over a variable number of time origins? 

    start_tau: The first lag time to calculate; subsequent lag times are calculated to be start_tau + i * delta_tau, where i = 1, 2... 
    and delta_tau is the step size of the lag time. 

    delta_tau: The step size of the lag times. By default, this is 1. For a delta_tau of 1, the lag times will be 1, 2, 3, ..., max_tau

    '''
    def __init__(self, energy_diff_array, max_tau = None, start_tau = 1, delta_tau = 1, block_size = 1):
        self.n_frames = energy_diff_array.shape[1] # Number of columns in the energy difference array (counting from start_t0)
        self.start_t0 = 0 # The first time origin: this is 0 since the array was only started from start_t0
        self.start_tau = start_tau # The first lag time 
        self.delta_tau = delta_tau # Step size in the lag times tau. 
        self.block_size = block_size # Block size in calculating the mean energy array 

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
        self.n_atoms = energy_diff_array.shape[0] # Number of atoms = number of rows in energy difference array  

        # Mean potential energy per atom
        if block_size<=1:
             mean_energy_single_col = np.mean(energy_diff_array, axis=1)
        else:
            mean_energy_single_col = np.mean(arr[:(len(arr)//block_size)*block_size].reshape(-1,block_size), axis=1)

        # Broadcast into a repeated array to match the shape of energy_diff_array
        mean_energy = np.transpose([mean_energy_single_col] * self.n_frames)
        # Array of energy fluctuations ( input for C(t) )
        self.energ_fluc = energy_diff_array - mean_energy

    def e_tau_e_t0(self, f_energ_t0, f_energ_t):
        ''' This calculates the product of the fluctuations in the energy:
        f_energ_t0 (fluctuation in the energy at a particular time origin)
        and f_energ_t (at a particular lag time), 
        averaged over all the particles in the system (say, n_atoms). 
        This returns the a single value.

        f_energ_t0: NumPy array of fluctuations in the energy of n_atoms at a particular time origin;
                        should have a shape of (n_atoms, dim), where dim is the number of dimensions
                        for which the TCF is being calculated.
        f_energ_t : NumPy array of fluctuations in the energy of n_atoms at time=t0+tau 
                        (where tau is a particular lag time);
                        should also have a shape of (n_atoms, dim)
        '''
        e_e0 = f_energ_t0 * f_energ_t

        # Get a single value (the numerator in the C(t) equation)
        return np.mean(e_e0, dtype=np.float64)


    def tcf_tau(self, tau):
        ''' This calculates the TCF (aka C(t), see 10.1103/PhysRevLett.93.023004 ), for a particular lag time tau, 
        averaged over all permitted time origins [ also averaged over all the particles in the system (say, n_atoms)]. 
        At this point, we already know the first time origin from which  
        which the TCF values will be calculated from (start_t0), and we know the number of time origins (n_origins)
        This returns the averaged product of the energy fluctuations, in the same units as the per atom energy values in
        the input data.

        tau: The lag time 

        Returns a single value 
        '''
        # TCF values of length n_origins for every time origin  
        tcf_origins = np.zeros(self.n_origins) 

        # Loop over all the time origins, starting from 
        # t0 = start_t0
        for idx in range(self.n_origins):
            t0 = idx + self.start_t0 # current time origin
            t = t0 + tau # time t = t0 + lag time 

            # Call the TCF function for t0 and tau
            tcf_origins[idx] = self.e_tau_e_t0(self.energ_fluc[:, t0], self.energ_fluc[:, t])

        # Get the average over all the permitted time origins
        return np.mean(tcf_origins, dtype=np.float64)

    def calculate_tcf(self):
        ''' Given a valid trajectory with velocities, calculate the TCF for various lag times. 
        The first lag time is given by start_tau, and subsequent values of the lag time
        are incremented in steps of delta_tau, up till the maximum lag time max_tau is reached. 

        The output numPy array will be of size (n_tau, 2), where n_tau corresponds to the number of lag times 
        calculated. The first column contains the lag times, with the other columns containing the TCF.
        '''
        # where we store the lag time and the associated TCF averaged over the time origins and atoms
        tcfList = [] # will be reshaped later 
        n_tau = 0 # Number of lag times 

        # Finding C(0) at t=0
        current_tau = 0;
        current_tcf = self.e_tau_e_t0(self.energ_fluc[:, 0], self.energ_fluc[:, 0])
        # Update the list of TCF values (1-D for now)
        tcfList.append(current_tau)
        tcfList.append(current_tcf)
        # Update the number of lag times 
        n_tau +=1; 

        # Loop over the lag times 
        for i_tau in range(self.start_tau, self.max_tau, self.delta_tau):
            n_tau += 1 
            current_tau = i_tau # Current lag time 
            # Get the TCF for this lag time over the desired time origins 
            current_tcf = self.tcf_tau(current_tau)
            # print("current tau: ", i_tau)
            # print("current vacf: ", current_vacf)
            # Update the list of TCF values (1-D for now)
            tcfList.append(current_tau)
            tcfList.append(current_tcf)

        # Get a numPy array of the lag times and the VACF values
        return np.array(tcfList).reshape( (n_tau, 2) )