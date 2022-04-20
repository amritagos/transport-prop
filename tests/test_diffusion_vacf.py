import pytest
import transportProp as trp
import numpy as np
import ase
from ase.atoms import Atom, Atoms # ASE stuff 
import pathlib
from . import aux

@pytest.fixture
def short_three_mols(): 
    ''' Reads in a LAMMPS trajectory of 7 steps, with three TIP4P/2005 water molecules, and returns
    a list of ASE Atoms objects with the coordinates and velocities. The time difference between each frame in the 
    trajectory is 1 ps (delta_t). The coordinates are in Angstroms, and are *wrapped*. The velocities of 
    every atom are provided. 
    Please note: this is actually a very short trajectory from the perspective of diffusion coefficient 
    calculations, using the Green-Kubo relation. Since every molecule consists of three atoms, 
    the ASE Atoms list returned actually contains the position of the O atoms only, and the *velocity*
    of the every corresponding molecule's *center of mass*. 
    This returns a VelocityAutoCorrelation object, initialized with the trajectory.
    '''
    p = aux.getpath('input/vel.lammpstrj') 
    traj = ase.io.read(p, format='lammps-dump-text', index=':') # Read in every step of the trajectory
    wat_traj = trp.misc.velocity_com_water_traj(traj)
    vacf_three_mol =  trp.VelocityAutoCorrelation(wat_traj)
    return vacf_three_mol

def test_vacf_all_atoms_t_t0(short_three_mols):
    ''' Tests the value of the VACF (for velocities in the x, y and z dimensions), 
    calculated for a particular time origin and lag time, averaged over all the particles
    (in this case, 3).  
    Takes in a VelocityAutoCorrelation object instantiated with the three particle trajectory 
    '''
    wat_traj = short_three_mols.traj # ASE Atoms list with velocities vx, vy, vz 
    # Given a time origin (the first frame, index 0 in wat_traj)
    # and a lag time of 1 delta_t (difference between time frames)
    vel_t0 = wat_traj[0].get_velocities() # Velocities of all atoms in the first frame
    vel_t = wat_traj[1].get_velocities() # Velocities of all atoms in the second frame
    vacf_t0_t = short_three_mols.v_tau_v_t0(vel_t0, vel_t) # VACF for one t0 and lag time, averaged over particles 
    # Check: 
    # should be  [ 4.37403330e-09 -1.50198503e-09 -2.53979207e-10]
    vacf_tt0_ref = np.array([ 4.37403330e-09, -1.50198503e-09, -2.53979207e-10])
    assert vacf_t0_t == pytest.approx(vacf_tt0_ref, 1e-6)
