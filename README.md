# Transport Property Calculations

This repository is meant to enable calculations of transport properties from equilibrium dynamics simulations. Presently, the source code is mostly in `Python`.

## Installation

We assume the existence of a reasonable build tool-chain. This means on an HPC you might need something like (on `elja`) `ml load GCC`.

### With `micromamba` or `conda` 

To install into a `micromamba` (or `conda`) environment, please follow the following steps. The `Python` dependencies are inside `environment.yml`. Create the environment using the following command, inside the top-level directory: 

```bash
micromamba create -f environment.yml
```
In order to activate the environment, use: 

```bash
micromamba activate trnsprtProp
```
where `trnsprtProp` is the name of the environment, specified in the `environment.yml` file. Henceforth, in order to activate the environment, just use the `micromamba activate trnsprtProp` command. If using `conda`, the equivalent commands are: 

```bash
conda env create -f environment.yml
conda activate trnsprtProp
```
To deactivate the environment, just use `micromamba deactivate`. 

Although this started out as a purely `Python` module (we started out using `flit`), we now use `pybind11`. Now we use `invoke` and `meson`. 

In order to install the code, run the following:

```bash
inv build # flit build
inv build --install
```

## Tests

To run tests (which are inside the `tests` directory), written with `pytest`, run the following command from the top-level directory: 

```bash
pytest
```

In order to debug tests using `pdb`, you can write the command `breakpoint()` inside the `Python` files (in `tests`) wherever you want to set a breakpoint. Then, run `pytest --pdb`. This will stop the code at the line where you put the `breakpoint()` command. 

To see more verbose output from `pytest`, including tests that pass, you can run `pytest -rA`.

## Running the code

This program can be run from the command-line using the following: 

```bash
transportProp
``` 

Help options can be accessed using:

```bash
transportProp --help 
```

Options for subcommands can be run, for instance: 

```bash
transportProp msd --help
```

## Validation 

Using the [VMD Diffusion Coefficient Tool](https://github.com/giorginolab/vmd_diffusion_coefficient), the MSD (mean-squared displacement) was calculated for a trajectory of mW water. In a `Jupyter Notebook`, the diffusion coefficient is calculated using linear regression. Please see the `validation` directory for more details. 
