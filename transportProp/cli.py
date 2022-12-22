from pathlib import Path 
import tomli # For reading the TOML file 
import numpy as np
from ase.atoms import Atoms # ASE stuff 
from ase.io import read, write, lammpsrun 

from . import structures 
from . import msd 
from . import vacf
from . import misc 
from . import tcf
from . import util
# from . import fastcpp

from typing import Optional
import typer 
from rich.console import Console

app = typer.Typer()

err_console = Console(stderr=True)

# Subcommands for MSD
@app.command()
def msd(
    tomlfile: Optional[Path] = typer.Argument(None, help="Path to the TOML file that contains options for the MSD."),
    traj: str = typer.Option(None, help="LAMMPS trajectory file, to be read in with ASE."),
    filetype: str = typer.Option(None, help="Trajectory file type, with the qualifier used by ASE, for instance, lammps-dump-text."),
    com: bool = typer.Option(True, "--com", help="When used, the centers of masses of the molecules are calculated. Currently hard-coded for water."),
    maxlag: int = typer.Option(None, min=1, help="Maximum lag time. When not specified, by default it is half the total number of steps."),
    firstorigin: int = typer.Option(None, min=0, help="A first time origin of 0 corresponds to the first configuration in the file. When not set by the user, this defaults to 0."),
    firstlag: int = typer.Option(None, min=1, help="The first lag time to calculate. When not specified, by default it is 1."),
    steplag: int = typer.Option(None, min=1, help="Step size in the lag times. When not specified, by default it is 1."),
    dim: str = typer.Option(None, help="This decides whether to calculate the MSD in the x, y, z or all dimensions. By default it calculates all dimensions. Allowed values are x, y, z, xyz or all."),
    outdir: Path = typer.Option(None, "-o", help="Path to the output directory. If not set, an output directory called output will be created.")
    ):
    """
    Calculates the mean-squared displacement. By default, the number of lag times is taken to be half of the 
    steps in the trajectory, thus ensuring that every MSD calculation averages over the same number of 
    time origins. 
    """
    if tomlfile is not None:
        if not tomlfile.exists():
            err_console.print("No TOML file.")
            raise typer.Abort()
        # Read in the TOML file
        msd_options = util.get_msd_data_options(tomlfile)
    else:
        msd_options = structures.MSDparams()

    # Optional arguments supercede the TOML file options
    if traj is not None:
          msd_options.trajectory = traj

    # Error handling for unspecified trajectory 
    if msd_options.trajectory is None:
        err_console.print("ERROR: You have not specified a trajectory file for the MSD calculation.")
        raise typer.Abort() # Abort 

    # Other CLI options for the MSD
    if filetype is not None:
        msd_options.trajectory_file_type = filetype
    if com:
        msd_options.use_center_of_mass = True
    if maxlag is not None:
        msd_options.max_lag_time = maxlag
    if firstorigin is not None:
        msd_options.first_time_origin = firstorigin
    if firstlag is not None:
        msd_options.first_lag_time = firstlag
    if steplag is not None:
        msd_options.step_size_lag_time = steplag
    if dim is not None:
        if dim!="x" and dim!="y" and dim!="z" and dim!="xyz" and dim!="all":
            err_console.print("Allowed values of dim are only: x, y, z, xyz, all.")
            raise typer.Abort()
        msd_options.dimension = dim
    if outdir is None:
        outdir = Path('output')

    # Calculate the MSD and write out the output files
    util.perform_msd_calc(msd_options, outdir)


# Subcommands for the TCF
@app.command()
def tcf(
    tomlfile: Optional[Path] = typer.Argument(None, help="Path to the TOML file that contains options for the TCF."),
    log0: str = typer.Option(None, help="LAMMPS log file, containing the energy gap for every timestep for the ground state."),
    log1: str = typer.Option(None, help="LAMMPS log file, containing the energy gap for every timestep for the excited state."),
    keystring: str = typer.Option(None, help="Qualifier or key string in the LAMMPS log file which denotes the energy gap value calculated for every timestep. For instance, c_solEnerg[0] could be a valid key string."),
    printdata: bool = typer.Option(True, "--print", help="When used, the timesteps and energy gap difference fluctuation arrays will be printed out from the log files."),
    timestring: str = typer.Option(None, help="Qualifier or key string in the LAMMPS log file for printing out the time steps."),
    maxlag: int = typer.Option(None, min=1, help="Maximum lag time. When not specified, by default it is half the total number of steps."),
    firstorigin: int = typer.Option(None, min=0, help="A first time origin of 0 corresponds to the first configuration in the file. When not set by the user, this defaults to 0."),
    firstlag: int = typer.Option(None, min=1, help="The first lag time to calculate. When not specified, by default it is 1."),
    steplag: int = typer.Option(None, min=1, help="Step size in the lag times. When not specified, by default it is 1."),
    outdir: Path = typer.Option(None, "-o", help="Path to the output directory. If not set, an output directory called output will be created.")
    ):
    """
    Calculates the solvation time correlation, which is based on the energy gap between the solute and solvent
    interactions in the ground and excited states. 
    """
    if tomlfile is not None:
        # Read in the TOML file
        tcf_options = util.get_tcf_data_options(tomlfile)
    else:
        tcf_options = structures.TCFparams()

    # print(fastcpp.add(3,4))

    # Optional arguments supercede the TOML file options
    # Ground state logfile
    if log0 is not None:
          tcf_options.log_ground_state = log0
    # Excited state logfile
    if log1 is not None:
          tcf_options.log_excited_state = log1
    # Key string for the computed energy gap values in the log files
    if keystring is not None:
        tcf_options.energy_gap_key_string = keystring

    # Error handling for unspecified ground state logfile 
    if tcf_options.log_ground_state is None:
        err_console.print("ERROR: You have not specified a log file (for the ground state), for the TCF calculation.")
        raise typer.Exit(code=1) # Exit with error  
    # Error handling for unspecified excited state logfile 
    if tcf_options.log_excited_state is None:
        err_console.print("ERROR: You have not specified a log file (for the excited state), for the TCF calculation.")
        raise typer.Exit(code=1) # Exit with error
    # Error handling for unspecified key string in the log file  
    if tcf_options.energy_gap_key_string is None:
        err_console.print("ERROR: You have not specified a key string for reading energy gap values, for the TCF calculation.")
        raise typer.Exit(code=1) # Exit with error

    # Other CLI options for the TCF
    if maxlag is not None:
        tcf_options.max_lag_time = maxlag
    if firstorigin is not None:
        tcf_options.first_time_origin = firstorigin
    if firstlag is not None:
        tcf_options.first_lag_time = firstlag
    if steplag is not None:
        tcf_options.step_size_lag_time = steplag
    if outdir is None:
        outdir = Path('output')
    if timestring is not None:
        tcf_options.timestep_key_string = timestring
    if not printdata:
        printdata = False

    # Calculate the TCF and write out the output files
    util.perform_tcf_calc(tcf_options, outdir, printdata)

if __name__ == "__main__":
    app()