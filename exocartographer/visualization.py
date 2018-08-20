# Import statements
import numpy as np
import matplotlib.pyplot as plt
from analytic_lc import analytic_lightcurve_2
import healpy as hp
from astropy.coordinates import Angle
import astropy.units as u

# Inputs for the light-curve
a = Angle((np.pi/2)*u.rad)
lat = Angle((np.pi/4)*u.rad)
colat_1 = a-lat
lng = Angle((0.8)*u.rad)

inclination = 0
obliquity = 0
sol_phase = 0

p_rotation = 23.934
p_orbit = 365.256363 * 24.0
w_rot = 2*np.pi/p_rotation
w_orb = 2*np.pi/p_orbit

nside=4
times = np.linspace(start=0.0, stop=24.0, num=1400)

analytic_1 = analytic_lightcurve_2(colat=colat_1, lng=lng, w_rot=w_rot, w_orb = w_orb,
                                   obq = obliquity, i=inclination, sol_phase=sol_phase,
                                   times=times, nside=nside)

# Plotting the light-curve
plt.plot(times, analytic_1)
plt.xlabel('time(h)')
plt.ylabel('reflectance')
plt.show()