import transportProp as trp
import numpy as np
import ase
from ase.atoms import Atom, Atoms # ASE stuff 
import pathlib
from ase.io import read, write, lammpsrun 

# # Read in the trajectory 
# wat_traj = read('MWwater-260', format='lammps-dump-text', index=':')

wat_traj = read('MWwater-260', format='lammps-dump-text', index=':')
# print(wat_traj[-1])

msd_mw = trp.MeanSquaredDisplacement(wat_traj, delta_tau = 10)
msdList = msd_mw.calculate_msd()

np.savetxt('msd-xyz.txt', msdList, delimiter=' ') 