# Declare exported functions
__all__ = [ 'getpath', 'get_slice', 'com_water_traj']

from pathlib import Path
from ase.build import molecule
from ase import Atoms, Atom

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