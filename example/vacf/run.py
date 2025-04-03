import transportProp as transportProp
import numpy as np
import ase
from ase.atoms import Atom, Atoms # ASE stuff 
import pathlib
from ase.io import read, write, lammpsrun 

from transportProp import structures, cli 

# # Read in the TOML file 
vacf_options = cli.get_vacf_data_options('input.toml')

# # Calculate the VACF and 
## write out the output files

cli.perform_vacf_calc(vacf_options)