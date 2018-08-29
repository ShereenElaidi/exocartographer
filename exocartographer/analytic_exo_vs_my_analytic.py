import numpy as np
from analytic_lc import analytic_lightcurve_2
from astropy.coordinates import Angle
import astropy.units as u
from exocartographer import viewgeom, kernel
import matplotlib.pyplot as plt


# Setting the input parameters
times = np.linspace(start=0.0, stop=10 * 24.0, num=14000)
p_rotation = 23.934
p_orbit = 365.256363 * 24.0
w_rot = 2 * np.pi/p_rotation
w_orb = 2 * np.pi/p_orbit
inclination = np.pi/2
obliquity = np.pi/6
sol_phase = np.pi/6
nside = 4
# Creating the point of the ``bright spot''
lat_1 = Angle((np.pi/4) * u.rad)
lng = Angle((-np.pi/2) * u.rad)

# Converting the latitude value into co-latitude values
a = Angle((np.pi/2) * u.rad)
colat_1 = a-lat_1

xi_0 = 0

trigvals = viewgeom(times=times, wrot=w_rot, worb=w_orb,
                    obq=obliquity, inc=inclination, xisol=sol_phase,
                    xi0=xi_0)
exocartographer = kernel(np.sin(colat_1), np.cos(colat_1),
                         np.sin(lng), np.cos(lng), trigA=trigvals)
analytic = analytic_lightcurve_2(colat=colat_1, lng=lng, w_rot=w_rot,
                                 w_orb=w_orb, obq=obliquity, i=inclination,
                                 sol_phase=sol_phase, times=times, nside=nside)

residuals = analytic - exocartographer
# Plotting the light-curves together over one day

f, (ax1, ax2) = plt.subplots(2,1, sharex=True, sharey=False)
ax1.plot(times, exocartographer, label='analytic - exocartographer')
ax1.plot(times, analytic, label='analytic - my code')
ax1.grid(color='#e6e8ed', linestyle='-', linewidth=2)
ax1.legend()
ax1.set_ylabel('reflectance')
ax1.set_xlabel('time(h)')

ax2.plot(times, residuals, label='residuals')
ax2.grid(color='#e6e8ed', linestyle='-', linewidth=2)
ax2.legend()
ax2.set_ylabel('residuals')
ax2.set_ylim(-0.20, 0.20)

f.tight_layout()
f.subplots_adjust(hspace=0)
plt.show()