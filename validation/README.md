# Sample Translational Diffusion Coefficient Calculation

The `mw-260-msd` directory contains text files containing the calculated MSD (mean-squared displacement) for an NVT simulation of 4096 mW water particles, at 260 K. 

## Simulation Details

 - The thermostat used was the Nosé-Hoover thermostat, and the simulation ran for almost 1 μs. Every frame written in the trajectory file (a text LAMMPS trajectory file) was spaced 1 ns apart. 
 - Time lags were taken from 1 to half the trajectory length, with step sizes of 10 (in ns), ensuring that for every time lag, time averaging is done over the same number of time origins. 
 - The coordinates in the trajectory file were *unwrapped*. 
 - The output MSD produced by the [VMD Diffusion Coefficient Tool](https://github.com/giorginolab/vmd_diffusion_coefficient). 
 - The time lags in the files are in nanoseconds. The MSD is in Å².  