# Declare exported functions
__all__ = [ 'get_msd_data_options', 'perform_msd_calc']

import pathlib 
import tomli # For reading the TOML file 
import numpy as np
from ase.atoms import Atoms # ASE stuff 
from ase.io import read, write, lammpsrun 

from . import structures 
from . import msd 
from . import vacf
from . import misc 

def get_msd_data_options(toml_filename):
    ''' Read in the options for performing a mean-squared displacement calculation.

    toml_filename: String with the name of the TOML file containing options for the MSD calculation.
    Returns a structures.MSDparams object whose members have fields with the MSD options 
    such as the trajectory name, first lag time value, etc. 
    '''
    p = pathlib.Path(toml_filename) # Path object for the TOML file

    # Read in the TOML file (as a binary file) 
    with p.open('rb') as f:
        data = tomli.load(f)

    # Get the options for [msd] in the TOML file
    msd_options = structures.MSDparams(**data.get('msd'))
    return msd_options


def perform_msd_calc(msd_options):
    '''
    This function takes in a structures.MSDparams object, containing options for performing the MSD
    calculation. 

    The trajectory is then actually read in, to obtain a list of ASE Atoms objects corresponding to each configuration. 
    Then we create the MeanSquaredDisplacement object and calculate the MSD. 

    Write out the output files relative to the current directory. 
    '''
    # If the file type has not been provided by the user, then the file format will be inferred 
    # by ASE from the filename
    if msd_options.trajectory_file_type is None:
        # Read in all the steps 
        atoms_traj = read(msd_options.trajectory, index=':')
    else:
        atoms_traj = read(msd_options.trajectory, format=msd_options.trajectory_file_type, index=':')

    # Make output directory
    out_dir = 'output'
    path = pathlib.Path(out_dir) # TODO: time stamp?
    path.mkdir(parents=True, exist_ok=True) 

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
        np.savetxt(out_dir+'/msd-'+'x'+'.txt', msdList, delimiter=' ') 
        ## Y dimension
        msd_obj = msd.MeanSquaredDisplacement(atoms_traj, max_tau = msd_options.max_lag_time, 
            start_t0 = msd_options.first_time_origin, start_tau = msd_options.first_lag_time, 
            delta_tau = msd_options.step_size_lag_time, dim = 'y')
        # Calculate the MSD
        msdList = msd_obj.calculate_msd()
        # Write out to file
        np.savetxt(out_dir+'/msd-'+'y'+'.txt', msdList, delimiter=' ') 
        ## Z dimension
        msd_obj = msd.MeanSquaredDisplacement(atoms_traj, max_tau = msd_options.max_lag_time, 
            start_t0 = msd_options.first_time_origin, start_tau = msd_options.first_lag_time, 
            delta_tau = msd_options.step_size_lag_time, dim = 'z')
        # Calculate the MSD
        msdList = msd_obj.calculate_msd()
        # Write out to file
        np.savetxt(out_dir+'/msd-'+'z'+'.txt', msdList, delimiter=' ') 
        ## X, Y, Z dimensions combined 
        msd_obj = msd.MeanSquaredDisplacement(atoms_traj, max_tau = msd_options.max_lag_time, 
            start_t0 = msd_options.first_time_origin, start_tau = msd_options.first_lag_time, 
            delta_tau = msd_options.step_size_lag_time, dim = 'xyz')
        # Calculate the MSD
        msdList = msd_obj.calculate_msd()
        # Write out to file
        np.savetxt(out_dir+'/msd-'+'xyz'+'.txt', msdList, delimiter=' ') 
    else:
        ## the dimension provided by the user 
        msd_obj = msd.MeanSquaredDisplacement(atoms_traj, max_tau = msd_options.max_lag_time, 
            start_t0 = msd_options.first_time_origin, start_tau = msd_options.first_lag_time, 
            delta_tau = msd_options.step_size_lag_time, dim = msd_options.dimension)
        # Calculate the MSD
        msdList = msd_obj.calculate_msd()
        # Write out to file
        np.savetxt(out_dir+'/msd-'+msd_options.dimension+'.txt', msdList, delimiter=' ') 


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