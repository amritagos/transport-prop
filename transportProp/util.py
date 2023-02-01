# Declare exported functions
__all__ = [ 'get_msd_data_options', 'perform_msd_calc', 'get_vacf_data_options', 'perform_vacf_calc',
'get_tcf_data_options', 'perform_tcf_calc']

import pathlib 
import tomli # For reading the TOML file 
import numpy as np
import math 
from ase.atoms import Atoms # ASE stuff 
from ase.io import read, write, lammpsrun 
import lammps_logfile 

import typer 
from rich.console import Console

from . import structures 
from . import msd 
from . import vacf
from . import misc 
from . import tcf
from . import fastcpp

err_console = Console(stderr=True)

def get_msd_data_options(toml_filename):
    ''' Read in the options for performing a mean-squared displacement calculation.

    toml_filename: Path object for TOML file containing options for the MSD calculation.
    Returns a structures.MSDparams object whose members have fields with the MSD options 
    such as the trajectory name, first lag time value, etc. 
    '''

    # Read in the TOML file (as a binary file) 
    with toml_filename.open('rb') as f:
        data = tomli.load(f)

    # Get the options for [msd] in the TOML file
    msd_options = structures.MSDparams(**data.get('msd'))
    return msd_options


def perform_msd_calc(msd_options, output_path):
    '''
    This function takes in a structures.MSDparams object, containing options for performing the MSD
    calculation. 

    The trajectory is then actually read in, to obtain a list of ASE Atoms objects corresponding to each configuration. 
    Then we create the MeanSquaredDisplacement object and calculate the MSD. 

    Write out the output files given the output Path object 
    '''
    # If the file type has not been provided by the user, then the file format will be inferred 
    # by ASE from the filename
    if msd_options.trajectory_file_type is None:
        # Read in all the steps 
        atoms_traj = read(msd_options.trajectory, index=':')
    else:
        atoms_traj = read(msd_options.trajectory, format=msd_options.trajectory_file_type, index=':')

    # Make output directory
    output_path.mkdir(parents=True, exist_ok=True) 

    # Get the center of masses 
    # TODO: make this more general
    # Hard-coded for water 
    if msd_options.use_center_of_mass:
        atoms_traj = misc.com_water_traj(atoms_traj)

    # For all files to be written out: 
    if msd_options.dimension == 'all':
        ## X dimension
        msd_obj = msd.MeanSquaredDisplacement(atoms_traj, max_tau = msd_options.max_lag_time, 
            start_t0 = msd_options.first_time_origin, start_tau = msd_options.first_lag_time, 
            delta_tau = msd_options.step_size_lag_time, dim = 'x')
        # Calculate the MSD
        msdList = msd_obj.calculate_msd()
        # Write out to file
        np.savetxt(str(output_path)+'/msd-'+'x'+'.txt', msdList, delimiter=' ') 
        ## Y dimension
        msd_obj = msd.MeanSquaredDisplacement(atoms_traj, max_tau = msd_options.max_lag_time, 
            start_t0 = msd_options.first_time_origin, start_tau = msd_options.first_lag_time, 
            delta_tau = msd_options.step_size_lag_time, dim = 'y')
        # Calculate the MSD
        msdList = msd_obj.calculate_msd()
        # Write out to file
        np.savetxt(str(output_path)+'/msd-'+'y'+'.txt', msdList, delimiter=' ') 
        ## Z dimension
        msd_obj = msd.MeanSquaredDisplacement(atoms_traj, max_tau = msd_options.max_lag_time, 
            start_t0 = msd_options.first_time_origin, start_tau = msd_options.first_lag_time, 
            delta_tau = msd_options.step_size_lag_time, dim = 'z')
        # Calculate the MSD
        msdList = msd_obj.calculate_msd()
        # Write out to file
        np.savetxt(str(output_path)+'/msd-'+'z'+'.txt', msdList, delimiter=' ') 
        ## X, Y, Z dimensions combined 
        msd_obj = msd.MeanSquaredDisplacement(atoms_traj, max_tau = msd_options.max_lag_time, 
            start_t0 = msd_options.first_time_origin, start_tau = msd_options.first_lag_time, 
            delta_tau = msd_options.step_size_lag_time, dim = 'xyz')
        # Calculate the MSD
        msdList = msd_obj.calculate_msd()
        # Write out to file
        np.savetxt(str(output_path)+'/msd-'+'xyz'+'.txt', msdList, delimiter=' ') 
    else:
        ## the dimension provided by the user 
        msd_obj = msd.MeanSquaredDisplacement(atoms_traj, max_tau = msd_options.max_lag_time, 
            start_t0 = msd_options.first_time_origin, start_tau = msd_options.first_lag_time, 
            delta_tau = msd_options.step_size_lag_time, dim = msd_options.dimension)
        # Calculate the MSD
        msdList = msd_obj.calculate_msd()
        # Write out to file
        np.savetxt(str(output_path)+'/msd-'+msd_options.dimension+'.txt', msdList, delimiter=' ') 


def get_vacf_data_options(toml_filename):
    ''' Read in the options for performing a velocity autocorrelation calculation.

    toml_filename: String with the name of the TOML file containing options for the VACF calculation.
    Returns a structures.VACFparams object whose members have fields with the VACF options 
    such as the trajectory name, first lag time value, etc. 
    '''
    p = pathlib.Path(toml_filename) # Path object for the TOML file

    # Read in the TOML file (as a binary file) 
    with p.open('rb') as f:
        data = tomli.load(f)

    # Get the options for [msd] in the TOML file
    vacf_options = structures.VACFparams(**data.get('vacf'))
    return vacf_options

def perform_vacf_calc(vacf_options):
    '''
    This function takes in a structures.VACFparams object, containing options for performing the VACF
    calculation. 

    The trajectory is then actually read in, to obtain a list of ASE Atoms objects corresponding to each configuration.
    The velocities must be present in the trajectory! 
    Then we create the VelocityAutoCorrelation object and calculate the VACF for vx*vx, vy*vy, vz*vz. 

    Write out the output files relative to the current directory. 
    '''
    # If the file type has not been provided by the user, then the file format will be inferred 
    # by ASE from the filename
    if vacf_options.trajectory_file_type is None:
        # Read in all the steps 
        atoms_traj = read(vacf_options.trajectory, index=':')
    else:
        atoms_traj = read(vacf_options.trajectory, format=vacf_options.trajectory_file_type, index=':')

    # Make output directory
    out_dir = 'output'
    path = pathlib.Path(out_dir) # TODO: time stamp?
    path.mkdir(parents=True, exist_ok=True) # Overwrite?

    # Get the center of masses 
    # TODO: make this more general
    # Hard-coded for water 
    #
    # Returns a list of ASE Atoms objects, wherein each molecule is replaced by the O atom
    # coordinates, and the velocity of the center of mass. Assumes that the water is in the order O H H 
    if vacf_options.velocity_center_of_mass:
        atoms_traj = misc.velocity_com_water_traj(atoms_traj)

    ## the dimension provided by the user 
    vacf_obj = vacf.VelocityAutoCorrelation(atoms_traj, max_tau = vacf_options.max_lag_time, 
        start_t0 = vacf_options.first_time_origin, start_tau = vacf_options.first_lag_time, 
        delta_tau = vacf_options.step_size_lag_time)
    # Calculate the VACF
    vacfList = vacf_obj.calculate_vacf()

    # Write out to file
    np.savetxt(out_dir+'/vacf'+'.txt', vacfList, delimiter=' ', header = 'tau vacf_xx vacf_yy vacf_zz') 

def get_tcf_data_options(toml_filename):
    ''' Read in the options for performing a solvation time correlation calculation ( C(t) ).

    toml_filename: Path object for TOML file containing options for the TCF calculation.
    Returns a structures.TCFparams object whose members have fields with the TCF options 
    such as the log file names, first lag time value, etc. 
    '''

    # Read in the TOML file (as a binary file) 
    with toml_filename.open('rb') as f:
        data = tomli.load(f)

    # Get the options for [tcf] in the TOML file
    tcf_options = structures.TCFparams(**data.get('tcf'))
    return tcf_options

def perform_tcf_calc(tcf_options, output_path, skip_every, printdata):
    '''
    This function takes in a structures.TCFparams object, containing options for performing the TCF
    calculation. 

    Two log files are processed. One logfile should have the energy gaps for one state, 
    while the other one has energy gaps (per timestep) for the other state (but the same configuration). 
    Currently there is only support for LAMMPS log files. The key string or compute name for which 
    the energy gap is saved in the log file is required.  
    Then we create the TimeCorrelation object and calculate the unnormalized TCF. 

    Write out the output files relative to the current directory if no Path is given. 
    '''
    # Read the log files for the ground and excited states
    log_excited = lammps_logfile.File(tcf_options.log_excited_state)
    log_ground = lammps_logfile.File(tcf_options.log_ground_state)

    # Get the numPy arrays (of length nsteps) for the energy gaps
    # for the ground and excited states 
    energ_ground = log_ground.get(tcf_options.energy_gap_key_string)
    energ_excited = log_excited.get(tcf_options.energy_gap_key_string)
    # Get the time from the logfile 
    timestep = log_ground.get(tcf_options.timestep_key_string)
    timestep2 = log_excited.get(tcf_options.timestep_key_string) # should be the same as timestep

    # Check that timestep and timestep2 are the same: 
    if not np.array_equal(timestep, timestep2):
        err_console.print("ERROR: The time series for the ground and excited states are not of the same length.")
        raise typer.Exit(code=1) # Exit with error

    if skip_every is not None:
        energ_ground = np.array(energ_ground[0::int(skip_every)])
        energ_excited = np.array(energ_excited[0::int(skip_every)])
        timestep = np.array(timestep[0::int(skip_every)])

    # Make output directory
    output_path.mkdir(parents=True, exist_ok=True)

    # Energy gap difference E(t) between the ground and excited states  
    delta_energy = energ_excited - energ_ground

    # Number of frames
    n_frames = delta_energy.shape[0] # Number of times for which E(t) is available

    # Average over all frames
    mean_energy = np.mean(delta_energy)
    # Fluctuations (of size n_frames) w.r.t. to the mean for this trajectory 
    energ_fluc = delta_energy - mean_energy

    # Save the timestep data and energy fluctuations 
    # to the output directory
    if printdata:
        np.savetxt(str(output_path)+"/energ_fluc.csv", energ_fluc, delimiter=",")
        np.savetxt(str(output_path)+"/time.csv", timestep, delimiter=",")

    # ## Solvation TCF calculation, given the potential energy difference array  
    # tcf_obj = tcf.SolvationTimeCorrelation(energy_diff_array, max_tau = tcf_options.max_lag_time, 
    #     start_tau = tcf_options.first_lag_time, delta_tau = tcf_options.step_size_lag_time)
    # # Calculate the TCF
    # tcfList = tcf_obj.calculate_tcf()

    # By default, max_tau is the rounded down integer value of half of the total number of frames. 
    # Then, for every tau, the TCF is averaged over the same number of times.
    if tcf_options.max_lag_time is None:
        max_tau = math.floor(0.5*n_frames)
    elif tcf_options.max_lag_time >  math.floor(0.5*n_frames):
        max_tau = math.floor(0.5*n_frames) # to ensure averaging over the same number of time origins
    else:
        max_tau = tcf_options.max_lag_time

    tau_val, tcf_val, mean_t01 = fastcpp.time_corr_function(energ_fluc, timestep, max_tau, 
        tcf_options.first_time_origin, tcf_options.first_lag_time,
        tcf_options.step_size_lag_time)

    # Write out to file  
    header_string = 'tau\tsolv_tcf\tC(0)=' + str(mean_t01)
    np.savetxt(str(output_path)+'/tcf'+'.txt', np.column_stack((tau_val, tcf_val)), delimiter=' ', header = header_string) 