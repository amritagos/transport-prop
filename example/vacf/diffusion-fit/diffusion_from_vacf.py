'''
Script for obtaining D(t) as a "running integral" from the VACF.
In the file containing the computed VACF, the first column corresponds to the time (lag time).
The second, third and fourth columns correspond to the VACF computed from the velocities in the X, Y, Z directions, respectively
(i.e., vacf_xx, vacf_yy, vacf_zz).

Here, we assume that you have sampled enough data points.
Integrating using SciPy's simpson rule directly seems to give worse results. 
The trapezoidal rule seems to work better than the Simpson rule.

Contributors:
Rohit Goswami: <rog32@hi.is>
Amrita Goswami: <amrita@hi.is>
Elvar Örn Jónsson: <eojons@gmail.com>, <elvarorn@hi.is>

'''
import numpy as np
import pathlib 
import pylab as pl
import math
from scipy import integrate
import quadpy 

import pdb

## -----
# Variables: TODO (run in from the command line, or work into CLI, toml etc.)

delta_t_frame = 4 # 4 fs in this case 

## -----

## FUNCTIONS

# `values` should be sorted
def get_closest(array, values):
    # make sure array is a numpy array
    array = np.array(array)

    # get insert positions
    idxs = np.searchsorted(array, values, side="left")
    
    # find indexes where previous index is closer
    prev_idx_is_less = ((idxs == len(array))|(np.fabs(values - array[np.maximum(idxs-1, 0)]) < np.fabs(values - array[np.minimum(idxs, len(array)-1)])))
    idxs[prev_idx_is_less] -= 1
    
    return idxs

def extrapolate_y(xpoints, ypoints, x_idx, x, tol):
    '''
    Extrapolate the value of y given xpoints, x_idx
    '''
    if math.isclose(xpoints[x_idx], x, rel_tol=tol):
        return ypoints[x_idx]
    else:
        # Otherwise return the linearly interpolated value 
        # Get the bracketing values of x
        if x-xpoints[x_idx]>0:
            xlo_idx = x_idx 
            xhi_idx = x_idx+1
        else:
            xlo_idx = x_idx+1
            xhi_idx = x_idx
        # Linearly interpolated value
        xlo = xpoints[xlo_idx] 
        xhi =  xpoints[xhi_idx]
        ylo = ypoints[xlo_idx]
        yhi = ypoints[xhi_idx]
        # y between (xlo, ylo) and (xhi, yhi)
        y = ylo + (x-xlo)*(yhi-ylo)/(xhi-xlo)
        return y

# Multiple values of x 
def gen_function(x, xpoints, ypoints, tol=1e-9):
    '''
    Several numerical integration techniques require a function. If the data is in the form of discrete samples,
    and assuming that the frequency of sampling is high enough, then any intermediate point in between two points 
    can be determined using linear interpolation.

    x : x values for which f(x) is desired 
    xpoints: discrete sample data x values. Assume that these are ordered, in ascending order. 
    ypoints: samples of f(x) available

    '''

    # Find the index of the nearest x point to the given x
    #pdb.set_trace()
    x_ids = get_closest(xpoints, x)
    
    y_out = x # Same shape etc. as x; now fill with new values
    len_x = np.shape(x)[0]
    
    for i in range(len_x):
        x_idx = x_ids[i]
        xVal = x[i]
        y_out[i] = extrapolate_y(xpoints, ypoints, x_idx, xVal, tol)
        
    return y_out

## --------------------------------

# Path to the VACF filename
filename = '../output/vacf.txt'
p = pathlib.Path(filename) # Path object for the VACF file 

data = np.loadtxt(p, skiprows=1)

time = data[:,0]*delta_t_frame # time in fs 
vacf_xx = data[:,1]
vacf_yy = data[:,2]
vacf_zz = data[:,3]

# Set up the output Diffusion coefficient arrays etc. 

n_points = np.shape(time)[0] # Number of sample points 

diff_t = np.zeros(( int(n_points-3),4)) # Integrate from the third point onwards 

# Get running average (integral)
for i in range(3, n_points):
    diff_t[i-3, 0] = time[i] # Current t
    # -----------------
    # SciPy's quad integrate (fails)
    # # Integrate upto t
    # integral_xx = integrate.quad(gen_function, 0, time[i], args=(vacf_xx, time))
    # integral_yy = integrate.quad(gen_function, 0, time[i], args=(vacf_yy, time))
    # integral_zz = integrate.quad(gen_function, 0, time[i], args=(vacf_zz, time))
    # ----
    # SIMPSONS (DOES NOT WORK WELL)
    # integral_xx = np.abs(integrate.simpson( vacf_xx[0:i], time[0:i] ))
    # integral_yy = np.abs(integrate.simpson( vacf_yy[0:i], time[0:i] ))
    # integral_zz = np.abs(integrate.simpson( vacf_zz[0:i], time[0:i] ))
    # ----
    # TRAPEZOIDAL (WORKS FINE)
    integral_xx = np.abs(integrate.trapezoid( vacf_xx[0:i], time[0:i] ))
    integral_yy = np.abs(integrate.trapezoid( vacf_yy[0:i], time[0:i] ))
    integral_zz = np.abs(integrate.trapezoid( vacf_zz[0:i], time[0:i] ))

    # Update the instantaneous diffusion coeff
    diff_t[i-3, 1] = integral_xx # D_xx {take first element in the case of quad}
    diff_t[i-3, 2] = integral_yy # D_yy
    diff_t[i-3, 3] = integral_zz # D_zz
    # # -----------
    # # Gauss Patterson scheme
    # scheme = quadpy.c1.gauss_patterson(5)
    # integral_xx = scheme.integrate( lambda x: gen_function(x, time, vacf_xx) , [0, time[i]] )
    # integral_yy = scheme.integrate( lambda x: gen_function(x, time, vacf_yy) , [0, time[i]] )
    # integral_zz = scheme.integrate( lambda x: gen_function(x, time, vacf_zz) , [0, time[i]] )

    # # Update the instantaneous diffusion coeff
    # diff_t[i-3, 1] = integral_xx # D_xx {take first element in the case of quad}
    # diff_t[i-3, 2] = integral_yy # D_yy
    # diff_t[i-3, 3] = integral_zz # D_zz
    # # -----------

# Make output directory
out_dir = 'output'
path = pathlib.Path(out_dir) # TODO: time stamp?
path.mkdir(parents=True, exist_ok=True) # Overwrite?
# Write out to file
np.savetxt(out_dir+'/diff'+'.txt', diff_t, delimiter=' ', header = 't D_xx D_yy D_zz') 

# for i in range(len(D)):
#     D[i] += (vacf_xx[:i+1].sum() + vacf_yy[:i+1].sum() + vacf_zz[:i+1].sum()) / 3.0
# pl.plot(D)
# pl.show()