# Declare exported functions
__all__ = [ 'getpath', 'get_slice', 'com_water_traj', 'v_com_molecule', 'velocity_com_water_traj']

from pathlib import Path
from ase.build import molecule
from ase import Atoms, Atom
import numpy as np
from scipy.optimize import curve_fit

def getpath(*a):
    ''' Get the directory in which the current file lives 
    '''
    # Package root
    d = Path(__file__).parent.resolve()
    return d.joinpath(*a)

def get_slice(a, dim):
    ''' Returns a sliced array, given an array a. 
    dim should be a string which is either 'x', 'y' or 'z'.
    '''
    if dim=='x':
        return a[: , 0]
    elif dim=='y':
        return a[: , 1]
    elif dim=='z':
        return a[:, 2]
    else:
        # for the option dim='xyz'
        return a 

def com_water_traj(atoms_traj):
    ''' Given a list of ASE Atoms objects, (read in from a trajectory)
    and assuming that it is for water and is in the order O H H , 
    get the center of mass of each molecule, and return a new Atoms object.

    We assume that there are no other atom types in the trajectory, and that the water
    molecules are perfectly ordered. 

    TODO: Figure out a way to make this more general 
    '''
    n_atoms = len(atoms_traj[0]) # Number of atoms in every configuration / frame (assume unchanged in time)
    wat = molecule('H2O') # Its positions will be set to the current molecule. 
    new_traj = [] # New trajectory containing the centers of masses of the water molecules 

    for atoms in atoms_traj:
        pos = atoms.get_positions() # positions of all O H H atoms in the current frame
        current_frame = Atoms() # center of masses of the current molecule 
        for i in range(0,n_atoms,3):
            wat.set_positions(pos[i:i+3, :])
            com = wat.get_center_of_mass() # Center of mass of the current molecule 
            current_frame += Atom('O', position = com) # Add com to the current frame 
        new_traj.append(current_frame) # Append frame to the trajectory 

    # Return the list of the Atoms objects with COM 
    return new_traj

def v_com_molecule(masses, v_atoms):
    '''
    Calculates the velocity of the center of mass of a molecule.

    masses: list of float; Atomic masses in atomic units. This must be in the same order as the velocity array
    
    v_atoms: numPy array of size (n_atoms, dim) where n_atoms is the number of atoms in the molecule
    and dim is the number of dimensions 

    Returns the a numPy array with velocity of the center of mass of a molecule, of size dim. 
    '''
    n_atoms = np.shape(v_atoms)[0]
    dim = np.shape(v_atoms)[1]
    vel_mol = np.zeros(dim)
    
    for iatom in range(n_atoms):
        vel_mol[:] += masses[iatom]*v_atoms[iatom, :] 
    
    return vel_mol/np.sum(masses)


def velocity_com_water_traj(atoms_traj):
    ''' Given a list of ASE Atoms objects, (read in from a trajectory)
    and assuming that it is for water and is in the order O H H , 
    get the velocity of the center of mass of each molecule, and return a new list of Atoms.

    We assume that there are no other atom types in the trajectory, and that the water
    molecules are perfectly ordered. 

    NOTE: The positions are not the COM of the molecules and are the positions of the O atoms only.
    This was done since the COM function does not take the minimum image convention into account.

    Additionally, the velocities are set in a rather circuitous way. This is because ASE Atoms objects store 
    momenta, not the velocities directly. When reading in LAMMPS trajectories, the symbols (and consequently the atomic masses)
    may not be correct. Therefore, it is important to obtain / set the velocities using the get_velocities() and set_velocities()
    functions, respectively. 

    TODO: Figure out a way to make this more general!! 
    '''
    n_atoms = len(atoms_traj[0])
    wat = molecule('H2O')
    wat_masses = wat.get_masses()
    new_traj = []
    v_frame = [] # velocities of molecules in one frame at a particular time 

    for atoms in atoms_traj:
        Oatoms = atoms[::3] 
        vel_atoms = atoms.get_velocities()
        for i in range(0, n_atoms, 3):
            v_com = v_com_molecule(wat_masses, vel_atoms[i:i+3,:])
            v_frame.append(v_com)
        # Handle the velocities
        v_frame = np.array(v_frame)
        v_frame = np.reshape(v_frame, (len(Oatoms),3))
        Oatoms.set_velocities(v_frame)
        new_traj.append(Oatoms)
        v_frame = []

    # Return the list of the Atoms objects with the velocities 
    return new_traj

def get_energy_difference(atoms_traj0, atoms_traj1, key_value, start_t0):
    ''' Given a list of ASE Atoms objects, (read in from two trajectories)
    and assuming the order of atoms is the same in each, corresponding to the same configuration,
    process the potential energies per atom and return a numPy array of differences in the potential
    energy per atom (every column corresponds to a particular time step).  

    We will skip all time steps before the first time origin. 
    key_value refers to the heading of the column for the per atom potential energies
    in the LAMMPS trajectory files. 

    TODO: Figure out a way to make this more general!! 
    '''
    n_atoms = len(atoms_traj0[0]) # assume constant
    total_steps = len(atoms_traj0) # total number of configurations in the trajectories
    energy0 = atoms_traj0[start_t0].arrays[key_value].ravel() # PE per atom in traj0 (n_atoms, 1)
    energy1 = atoms_traj1[start_t0].arrays[key_value].ravel() # PE per atom in traj1 (n_atoms, 1)

    for istep in range(start_t0+1, total_steps):
        # Add columns corresponding to each configuration at every time step 
        energy0 = np.column_stack([energy0, atoms_traj0[istep].arrays[key_value].ravel()])
        energy1 = np.column_stack([energy1, atoms_traj1[istep].arrays[key_value].ravel()])

    # Return the numPy array 
    return energy0-energy1

def quadratic(x, a, b, c):
    """The quadratic function to fit."""
    return a * x**2 + b * x + c

def two_part_quadratic_fit(x, y):
    """Perform a two-part quadratic fitting procedure on the given data."""
    
    # Find the maximum value and its index
    max_index = np.argmax(y)
    max_x = x[max_index]
    
    # Fit a quadratic using 10 points on either side of the maximum value
    left_index = max(0, max_index - 10)
    right_index = min(len(x), max_index + 11)
    quadratic_1_x = x[left_index:right_index]
    quadratic_1_y = y[left_index:right_index]
    popt_1, pcov_1 = curve_fit(quadratic, quadratic_1_x, quadratic_1_y)
    
    # Find the x value closest to the apex of the first quadratic fit
    apex_x = -popt_1[1] / (2 * popt_1[0])
    apex_index = np.argmin(np.abs(x - apex_x))
    
    # Fit a quadratic using 10 points on either side of the apex of the first quadratic fit
    left_index = max(0, apex_index - 10)
    right_index = min(len(x), apex_index + 11)
    quadratic_2_x = x[left_index:right_index]
    quadratic_2_y = y[left_index:right_index]
    popt_2, pcov_2 = curve_fit(quadratic, quadratic_2_x, quadratic_2_y)
    
    return popt_1, popt_2