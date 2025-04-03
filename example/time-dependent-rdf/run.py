import transportProp as transportProp
from pathlib import Path

from transportProp import structures, util

## For testing with gdb
# Run the following: 
# gdb -args python run.py 

input_traj = 'dump/dump-500000.lammpstrj'
typeO = 1
typeIon = 3
cutoff = 15.0 
binwidth = 0.1 

input_path = Path(input_traj)

peak_pos, peak_heights = util.get_rdf_peakData_time(input_path, typeO, typeIon, binwidth, cutoff)