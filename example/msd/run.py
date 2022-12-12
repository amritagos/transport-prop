import transportProp as trp
import numpy as np
import ase
from ase.atoms import Atom, Atoms # ASE stuff 
import pathlib
from ase.io import read, write, lammpsrun 

from transportProp import structures, cli 

import typer 

app = typer.Typer()

# Subcommands for MSD
@app.command()
def msd(tomlfile: str):
    # Read in the TOML file
    msd_options = cli.get_msd_data_options(tomlfile)
    # Calculate the MSD and write out the output files
    cli.perform_msd_calc(msd_options)

if __name__ == "__main__":
    app()

# # # Read in the TOML file 
# msd_options = cli.get_msd_data_options('input.toml')

# # # Calculate the MSD and 
# ## write out the output files

# cli.perform_msd_calc(msd_options)