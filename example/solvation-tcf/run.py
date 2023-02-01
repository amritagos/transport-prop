import transportProp as trp
from pathlib import Path

from transportProp import structures, util

## For testing with gdb
# Run the following: 
# gdb -args python run.py 

# # Read in the TOML file 
tcf_options = util.get_tcf_data_options(Path('input.toml'))

# Do the calculation 
util.perform_tcf_calc(tcf_options, Path('output'), None, True)