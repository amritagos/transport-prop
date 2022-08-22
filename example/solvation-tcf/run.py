import transportProp as trp
import numpy as np
import ase
from ase.atoms import Atom, Atoms # ASE stuff 
import pathlib
from ase.io import read, write, lammpsrun 

from transportProp import structures, cli 

# # Read in the TOML file 
tcf_options = cli.get_tcf_data_options('input.toml')

# Do the calculation 
cli.perform_tcf_calc(tcf_options)