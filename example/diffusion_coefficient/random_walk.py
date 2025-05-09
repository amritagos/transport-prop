from pathlib import Path
import numpy as np
from numba import njit

@njit
def random_walk(n_steps, dz):
    for i in range(1, n_steps):
        # Choose step: +dz or -dz
        step = dz if np.random.rand() < 0.5 else -dz
        z_positions[i] = z_positions[i - 1] + step
    return z_positions

# Output folder 
output_dir = input_file_path = Path(__file__).parent / "output"
# Create the directory if it doesn't exist 
output_dir.mkdir(parents=True, exist_ok=True)

# Write trajectory to a LAMMPS dump-style file
output_filename = output_dir / "oxygen_random_walk.dump"

# One-dimensional random walk in the z-dimension
# Parameters
dt = 1.0  # time step in femtoseconds
dz = 2.0  # displacement in Angstrom per step (up or down)
n_steps = 5000  # total number of steps
atom_id = 1  # single atom id
atom_type = 1  # LAMMPS atom type

# Analytical diffusion coefficient
D_analytical = (dz**2) / (2 * dt)
# print("Analytical diffusion coefficient: {:.3f} Å²/fs".format(D_analytical))

# Initialize the trajectory
time = np.arange(0, n_steps * dt, dt)
z_positions = np.zeros(n_steps)
# x and y remain zero
x_positions = np.zeros(n_steps)
y_positions = np.zeros(n_steps)

# Seed random number generator for reproducibility
np.random.seed(42)

# Random walk simulation
random_walk(n_steps, dz)

with open(output_filename, "w") as f:
    for step in range(n_steps):
        f.write("ITEM: TIMESTEP\n")
        f.write(f"{step}\n")
        f.write("ITEM: NUMBER OF ATOMS\n")
        f.write("1\n")
        # Define box bounds arbitrarily (we only care about z)
        # Here, we set x and y from -10 to 10, and z bounds based on min and max of z_positions seen so far.
        zmin = np.min(z_positions[: step + 1]) - 1
        zmax = np.max(z_positions[: step + 1]) + 1
        f.write("ITEM: BOX BOUNDS pp pp pp\n")
        f.write("-10 10\n")
        f.write("-10 10\n")
        f.write(f"{zmin:.3f} {zmax:.3f}\n")
        f.write("ITEM: ATOMS id type x y z\n")
        f.write(
            f"{atom_id} {atom_type} {x_positions[step]:.3f} {y_positions[step]:.3f} {z_positions[step]:.3f}\n"
        )
print(f"Trajectory written to {output_filename}")

# # -------------------------------------
from ase.io import read
import transportProp as trp
from scipy.optimize import curve_fit
import matplotlib.pyplot as plt

max_tau = 400

# Read in the trajectory and process it
traj = read(output_filename, index=":")

msd_obj = trp.MeanSquaredDisplacement(
    traj, max_tau=max_tau, start_t0=0, start_tau=1, delta_tau=1, dim="z"
)
# Calculate the MSD
msdList = msd_obj.calculate_msd()
tau_val = msdList[:,0]
msd = msdList[:,1] 
# Write out to file
output_msd = output_dir / "msd_z.txt"
header = 'time_lag msd' # time lag and MSD value
np.savetxt(output_msd, msdList, delimiter=' ', header=header)

# ---------------------------------------
# Diffusion coefficient from a fit to the MSD 

# --- Fitting MSD = b + m * tau^n ---
# To avoid the short-time regime where the Einstein relation might not yet hold,
# set a lower time limit (t_min).
t_min = 100  # in fs, adjust this value as necessary for your system
t_max = max_tau
fit_mask = (tau_val >= t_min) & (tau_val < t_max)

# Prepare data for fitting
tau_fit = tau_val[fit_mask]
msd_fit = msd[fit_mask]

# Define the model function
def msd_model(tau, b, m, n):
    b = 0.0
    return b + m * tau**n

def linear_model(tau, b, m):
    return b + m*tau

# Provide an initial guess for the parameters: b, m, n
# For an ideal random walk, we expect b ~ 0, m ~ 4 (since MSD = 4 fs⁻¹*fs for our 1D case),
# and n ~ 1.
p0 = [0.0, 4.0, 1.0]
# p0 = [0.0, 4.0]

# Perform the non-linear fit
popt, pcov = curve_fit(msd_model, tau_fit, msd_fit, p0=p0)
b_fit, m_fit, n_fit = popt
# b_fit, m_fit = popt

# Calculate the fitted diffusion coefficient
# In 1D diffusion, the Einstein relation gives MSD = 2*D*tau.
# If the fit is good and n_fit ~ 1, then we identify m_fit = 2*D.
dim = 1 # one-dimensional case 
D_fit = m_fit / (2.0*dim)

print("Fitting MSD = b + m*tau^n")
print("Fitted parameters:")
print("  b      = {:.3f}".format(b_fit))
print("  m      = {:.3f}".format(m_fit))
# print("  n      = {:.3f}".format(n_fit))
print("Fitted diffusion coefficient D_fit = {:.3f} Å²/fs".format(D_fit))
print(f"{D_analytical=}")

# Plot the data and the fit for visualization
plt.figure(figsize=(8, 6))
plt.loglog(tau_fit, msd_fit, 'o', markersize=4, label=f'Data (t ≥ {t_min} fs)')
tau_plot = np.linspace(t_min, np.max(tau_fit), 200)
plt.loglog(tau_plot, msd_model(tau_plot, *popt), '-', label=f'Fit: n={n_fit}')
# plt.axhline(1.0)
# plt.loglog(tau_fit, np.gradient(msd_fit))
# plt.loglog(tau_plot, linear_model(tau_plot, *popt), '-', label=f'Fit: n=1')
plt.xlabel("Time (fs)")
plt.ylabel("MSD (Å²)")
plt.title("Fit of MSD vs. Time")
plt.legend()
plt.grid(True, which="both", ls="--")
plt.savefig(output_dir / "msd_fit_nonlinear.png")
plt.show()
