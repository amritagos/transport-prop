import pytest
import transportProp as trp
import numpy as np
import ase
from ase.atoms import Atom, Atoms # ASE stuff 
import pathlib
from . import aux

@pytest.fixture
def small_system(): 
    ''' Reads in two LAMMPS trajectories of 5 steps, both of which contain two TIP4P/2005 
    water molecules and one 'ground-state' ion, and returns
    a numPy array of the energy differences per atom, with columns corresponding to every 
    timestep and configuration 
    The time difference between each frame in the 
    trajectory is 10 ps (delta_t). The coordinates are in Angstroms, and are *wrapped*. 

    Please note: this is actually a very short trajectory from the perspective of the solvation correlation function
    This returns an ASE Atoms object, which must be processed further. The solvation correlation function requires
    the fluctuations in the energy differences as input.
    '''
    # State 0 ("Ground" state)
    p0 = aux.getpath('input/pe_energy_state0.lammpstrj') 
    traj0 = ase.io.read(p0, format='lammps-dump-text', index=':') # Read in every step of the trajectory
    # State 1 ("Excited" state)
    p1 = aux.getpath('input/pe_energy_state1.lammpstrj') 
    traj1 = ase.io.read(p1, format='lammps-dump-text', index=':') # Read in every step of the trajectory
    energy_diff = trp.misc.get_energy_difference(traj0, traj1, key_value='c_peratom', start_t0=0)
    return energy_diff

def test_energy_diff_array(small_system):
    ''' Tests that the energy difference has been computed correctly, for a system with one ion 
    and two water molecules. 
    Takes in a numPy array 
    '''
    diff_array = small_system # numPy array of differences 
    # Should have the shape (7,5)
    assert diff_array.shape == (7,5) 
    # # Given a time origin (the first frame, index 0 in wat_traj)
    # # and a lag time of 1 delta_t (difference between time frames)
    # vel_t0 = wat_traj[0].get_velocities() # Velocities of all atoms in the first frame
    # vel_t = wat_traj[1].get_velocities() # Velocities of all atoms in the second frame
    # vacf_t0_t = short_three_mols.v_tau_v_t0(vel_t0, vel_t) # VACF for one t0 and lag time, averaged over particles 
    # # Check: 
    # # should be  [ 4.37403330e-09 -1.50198503e-09 -2.53979207e-10]
    # vacf_tt0_ref = np.array([ 4.37403330e-09, -1.50198503e-09, -2.53979207e-10])
    # assert vacf_t0_t == pytest.approx(vacf_tt0_ref, 1e-6)
