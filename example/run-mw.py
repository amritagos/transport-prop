import transportProp as trp
import numpy as np
import ase
from ase.atoms import Atom, Atoms # ASE stuff 
import pathlib
from ase.io import read, write, lammpsrun 

from transportProp import structures, cli 

# # Read in the TOML file 
msd_options = cli.get_msd_data_options('input.toml')

# # Calculate the MSD and 
## write out the output files

cli.perform_msd_calc(msd_options)