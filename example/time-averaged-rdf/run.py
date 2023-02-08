import transportProp as trp
from pathlib import Path
import numpy as np

from transportProp import structures, util

## For testing with gdb
# Run the following: 
# gdb -args python run.py 

input_traj = 'dump/Fe2-equil.lammpstrj'
typeO = 1
typeIon = 3
cutoff = 15.0 
binwidth = 0.05 

input_path = Path(input_traj)

r_grid, rdf_avg = util.get_time_averaged_rdf(input_path, typeO, typeIon, binwidth, cutoff)

# Save to a file in an output folder 
header_string = 'r\tg(r)'
np.savetxt('output/rdf.txt', np.column_stack((r_grid, rdf_avg)), delimiter=' ', header = header_string)