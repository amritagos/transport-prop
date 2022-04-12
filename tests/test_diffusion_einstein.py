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

def test_traj(get_short_single_particle_traj):
    wat_traj = get_short_single_particle_traj # ASE Atoms list with unwrapped coordinates
    assert len(wat_traj) == 11
