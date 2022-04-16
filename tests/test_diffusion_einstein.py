import pytest
import transportProp as trp
import numpy as np
import ase
from ase.atoms import Atom, Atoms # ASE stuff 
import pathlib
from . import aux


@pytest.fixture
def get_short_single_particle_traj(): 
    ''' Reads in an XYZ trajectory of 11 steps, with a single mW (monatomic water) "particle", and returns
    a list of ASE Atoms objects with the coordinates. The time difference between each frame in the 
    trajectory is 1 ns (delta_t). The coordinates are in Angstroms, and are *unwrapped*. 
    Please note: this is actually a very short trajectory from the perspective of diffusion coefficient 
    calculations, using the Einstein relation (which involves doing a linear regression fit of 
    MSD data at various 'lag times'). For molecules, it would be preferable to use the coordinates 
    of the centers of masses.    
    '''
    p = aux.getpath('input/single-particle-short.xyz') 
    wat_traj = ase.io.read(p, index=':') # Read in every step of the trajectory
    return wat_traj

@pytest.fixture
def get_short_three_particle_traj(): 
    ''' Reads in an XYZ trajectory of 11 steps, with three mW (monatomic water) "particles", and returns
    a list of ASE Atoms objects with the coordinates. The time difference between each frame in the 
    trajectory is 1 ns (delta_t). The coordinates are in Angstroms, and are *unwrapped*. 
    Please note: this is actually a very short trajectory from the perspective of diffusion coefficient 
    calculations, using the Einstein relation (which involves doing a linear regression fit of 
    MSD data at various 'lag times'). For molecules, it would be preferable to use the coordinates 
    of the centers of masses.    
    '''
    p = aux.getpath('input/three-particles-short.xyz') 
    wat_traj = ase.io.read(p, index=':') # Read in every step of the trajectory
    return wat_traj

def test_msd_single_iatom_t_t0(get_short_single_particle_traj):
    ''' Tests the squared distance function between a given time origin and
    lag time for a single particle 
    '''
    wat_traj = get_short_single_particle_traj # ASE Atoms list with unwrapped coordinates
    # Given a time origin (the first frame, index 0 in wat_traj)
    # and a lag time of 1 delta_t (difference between time frames)
    # Calculate the squared distance for a single particle. 
    pos_t0 = wat_traj[0].get_positions() # Positions of all atoms in the first frame
    # Although here there is only one atom 
    pos_t = wat_traj[1].get_positions() # Positions of all atoms in the second frame
    # Get the position of the atom with index 0 at time origin t0  
    pos_t0_iatom = pos_t0[0 , :]
    # Get the positions of the atom with index 0 at time origin t
    pos_t_iatom = pos_t[0 , :]
    r2 = trp.sq_disp_iatom(pos_t0_iatom, pos_t_iatom) # Squared distance 
    # Check: 
    # should be 1482.2965976770176 A^2
    r2_ref = 1482.2965976770176
    assert r2 == pytest.approx(r2_ref, 1e-6)

def test_msd_all_atoms_t_t0(get_short_three_particle_traj):
    ''' Tests the value of the mean-squared displacement, 
    calculated for a particular time origin and lag time, averaged over all the particles
    (in this case, 3).  
    '''
    wat_traj = get_short_three_particle_traj # ASE Atoms list with unwrapped coordinates
    # Given a time origin (the first frame, index 0 in wat_traj)
    # and a lag time of 1 delta_t (difference between time frames)
    pos_t0 = wat_traj[0].get_positions() # Positions of all atoms in the first frame
    pos_t = wat_traj[1].get_positions() # Positions of all atoms in the second frame
    msd_t0_t = trp.tau_t0(pos_t0, pos_t) # MSD for one t0 and lag time, averaged over particles 
    # Check: 
    # should be 3109.92125715317 A^2
    msd_tt0_ref = 3109.92125715317
    assert msd_t0_t == pytest.approx(msd_tt0_ref, 1e-6)