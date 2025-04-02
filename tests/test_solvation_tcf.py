import pytest
import transportProp as trp
import numpy as np
import ase
from ase.atoms import Atom, Atoms # ASE stuff 
import pathlib
from . import aux

# @pytest.fixture
# def small_system_energy_diff_array(): 
#     ''' Reads in two LAMMPS trajectories of 5 steps, both of which contain two TIP4P/2005 
#     water molecules and one 'ground-state' ion, and returns
#     a numPy array of the energy differences per atom, with columns corresponding to every 
#     timestep and configuration 
#     The time difference between each frame in the 
#     trajectory is 10 ps (delta_t). The coordinates are in Angstroms, and are *wrapped*. 

#     Please note: this is actually a very short trajectory from the perspective of the solvation correlation function
#     This returns an ASE Atoms object, which must be processed further. The solvation correlation function requires
#     the fluctuations in the energy differences as input.
#     '''
#     # State 0 ("Ground" state)
#     p0 = aux.getpath('input/pe_energy_state0.lammpstrj') 
#     traj0 = ase.io.read(p0, format='lammps-dump-text', index=':') # Read in every step of the trajectory
#     # State 1 ("Excited" state)
#     p1 = aux.getpath('input/pe_energy_state1.lammpstrj') 
#     traj1 = ase.io.read(p1, format='lammps-dump-text', index=':') # Read in every step of the trajectory
#     energy_diff = trp.misc.get_energy_difference(traj0, traj1, key_value='c_peratom', start_t0=0)
#     return energy_diff

# @pytest.fixture
# def small_system(small_system_energy_diff_array): 
#     ''' Given the energy differences between two LAMMPS trajectories of 5 steps, both of which contain two TIP4P/2005 
#     water molecules and one 'ground-state' ion, returns a SolvationTimeCorrelation object.
#     The time difference between each frame in the 
#     trajectory is 10 ps (delta_t). The coordinates are in Angstroms, and are *wrapped*. 

#     Please note: this is actually a very short trajectory from the perspective of the solvation correlation function
#     This returns an ASE Atoms object, which must be processed further. The solvation correlation function requires
#     the fluctuations in the energy differences as input.
#     '''
#     tcf_small_sys = trp.SolvationTimeCorrelation(small_system_energy_diff_array)
#     return tcf_small_sys

# def test_energy_diff_array(small_system_energy_diff_array):
#     ''' Tests that the energy difference has been computed correctly, for a system with one ion 
#     and two water molecules. 
#     Takes in a numPy array 
#     '''
#     diff_array = small_system_energy_diff_array # numPy array of differences 
#     # Should have the shape (7,5)
#     assert diff_array.shape == (7,5) 

# def  test_tcf_all_atoms_t_t0(small_system):
#     ''' Tests the value of the unnormalized TCF , calculated for a particular time origin and lag time, 
#     averaged over all the particles (in this case, 2 water molecules and 1 ion).  
#     Takes in a SolvationTimeCorrelation object instantiated with the energy difference array.
#     '''
#     energy_fluc_array = small_system.energ_fluc # numPy array of energy difference fluctuations (n_atoms, n_frames)
#     # TCF for one t0 and lag time, averaged over particles
#     tcf_t0_t = small_system.e_tau_e_t0(energy_fluc_array[:,0], energy_fluc_array[:,1])  
#     # Check: 
#     tcf_tt0_ref = 11.926569490688 # should be  11.926569490688
#     assert tcf_t0_t == pytest.approx(tcf_tt0_ref, 1e-6)
