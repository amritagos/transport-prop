# Transport Property Calculations

This repository is meant to enable calculations of transport properties from equilibrium dynamics simulations. Presently, the source code is mostly in `Python`.

## Installation

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

To work with the `Python` module (since it is Python only), we have used `flit` to handle
builds. This supports editable
installs, which means that the package does not need to be reinstalled
explicitly on changes. 

So, in order to install the `Python` code, all you have to do is run the following command in the top-level directory:

```bash
pip install -e .
```

### With `venv`

If you would prefer to use a virtual environment using `venv`, run the following to create and activate a virtual environment: 

```bash
python3 -m venv envName
source envName/bin/activate
```

Install the dependencies using the `requirements.txt` file, using `python3 -m pip install -r requirements.txt` (check out the [docs](https://packaging.python.org/en/latest/guides/installing-using-pip-and-virtual-environments/) for Windows). 

To install the module, run the command `pip install -e .` in editable mode. 

## Tests

To run tests (which are inside the `tests` directory), written with `pytest`, run the following command from the top-level directory: 

```bash
pytest
```

In order to debug tests using `pdb`, you can write the command `breakpoint()` inside the `Python` files (in `tests`) wherever you want to set a breakpoint. Then, run `pytest --pdb`. This will stop the code at the line where you put the `breakpoint()` command. 

To see more verbose output from `pytest`, including tests that pass, you can run `pytest -rA`. 

## Validation 

Using the [VMD Diffusion Coefficient Tool](https://github.com/giorginolab/vmd_diffusion_coefficient), the MSD (mean-squared displacement) was calculated for a trajectory of mW water. In a `Jupyter Notebook`, the diffusion coefficient is calculated using linear regression. Please see the `validation` directory for more details. 