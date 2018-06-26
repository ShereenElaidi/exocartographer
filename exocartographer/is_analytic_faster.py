import time
from analytic_delta import analytic_delta
from numerically_integrate import compute_flux
import numpy as np

# SIMULATING AN ENTIRE PLANET'S LIGHTCURVE

# Array of Lat (theta) and Lng (phi)
theta = [0,np.pi/6, np.pi/5, np.pi/4, np.pi/3, np.pi/2,
         2*np.pi/3, 3*np.pi/4, 5*np.pi/6, np.pi, 7*np.pi/6,
         5*np.pi/4,4*np.pi/3, 5*np.pi/3, 7*np.pi/4, 11*np.pi/6]
phi = [0,np.pi/6, np.pi/5, np.pi/4, np.pi/3, np.pi/2,
         2*np.pi/3, 3*np.pi/4, 5*np.pi/6, np.pi, 7*np.pi/6,
         5*np.pi/4,4*np.pi/3, 5*np.pi/3, 7*np.pi/4, 11*np.pi/6]
print(theta)

# CALCULATING THE TIME DIFFERENCE
t0 = time.clock()
# Put code for analytic in here
t1 = time.clock()
analytic_time = t1-t0
print(analytic_time)

t2 = time.clock()
# Put code for numerical integration here
t3 = time.clock()
numerical_time = t3 - t2
print(numerical_time)

# CALCULATING THE ACCURACY OF THE ANALYTIC MAP
