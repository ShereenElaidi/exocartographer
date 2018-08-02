import time as t
import numpy as np
import matplotlib.pyplot as plt
from analytic_lc import analytic_lightcurve
from analytic_lc_2 import analytic_lc_2
from math import sqrt
import healpy as hp
from exocartographer.gp_map import draw_map
from exocartographer import IlluminationMapPosterior
from exocartographer.util import logit, inv_logit
from astropy_healpix import HEALPix
import pdb


# ==============DARK PLANET EVERYWHERE, ONE BRIGHT POINT==========

# -------Planet parameters----------------------------------------


# the unit here will be hours
p_rotation = 23.934
p_orbit = 365.256363 * 24.0

# the unit here will be rotations
# p_rotation = 1
# p_orbit = 365.256363

# -------Parameters-----------------------------------------------
time_data = np.arange(0.0, 10.0, 0.01)
# w_rot = 2*np.pi
# w_orb = (2*np.pi)/365.256363

w_rot = 2*np.pi/p_rotation
w_orb = 2*np.pi/p_orbit

obliq = 90*np.pi/180.0
i = np.pi/2
sol_phase = 1

# Bright point. Drawing the albedo map using HEALPIX. First,
# Gaussian process properties & resolution:

nside = 4
whitenoise_relative_amp = 0.2
length_scale = 30*np.pi/2
albedo_mean = 0.5
albedo_std = 0.2

while True:
    one_point_map = draw_map(nside, albedo_mean, albedo_std,
                             whitenoise_relative_amp, length_scale)
    if min(one_point_map) > 0 and max(one_point_map) <1:
        break

# Simulating a map that is all black everywhere except for
# one bright spot
bright_spot = 75
for i in range(0, one_point_map.size):
    if i == bright_spot:
        one_point_map[i] = 1.00
    else:
        one_point_map[i] = 0.00

# Viewing the map
# hp.mollview(one_point_map, title='One bright point', cmap='gist_gray')

# Obtaining the latitude and longitude of this bright point
map = HEALPix(nside=nside, order='nested')
pt_lng, pt_lat = map.healpix_to_lonlat([bright_spot])


# Assigning these points to the input parameters
lat = pt_lat
lng = pt_lng

# =============NUMERICAL LIGHTCURVE=============================

# Setting the orbital properties

# the unit here will be hours
p_rotation = 23.934
p_orbit = 365.256363 * 24.0

# the unit here will be rotations
# p_rotation = 1
# p_orbit = 365.256363

phi_orb = 3*np.pi/2 + sol_phase
inclination = np.pi/2
obliquity =0
phi_rot = 3*np.pi/2 + sol_phase

# Observation schedule
start_day = 50
cadence = p_rotation/4.
epoch_duration = p_orbit       # p_orbit for year, p_rot for 1 day
epoch_duration = p_rotation

times = np.linspace(start=50.0*24, stop=51.0*24.0, num=1400)
# Times:

# Measurement uncertainties
measurement_std = 0.001



# Inputs for the code:
lat, lng = map.healpix_to_lonlat([bright_spot])
w_rot = 2*np.pi/p_rotation
w_orb = 2*np.pi/p_orbit
inclination = np.pi/2
obliquity = 0
sol_phase = 1
p_rotation = 2*np.pi/w_rot
p_orb = 2*np.pi/w_orb

# DO NOT CHANGE THIS.
phi_orb = abs(5*(np.pi/3) + sol_phase)
phi_rot = abs(2*np.pi - sol_phase)


# Obtaining the numerical lightcurve
truth = IlluminationMapPosterior(times, np.zeros_like(times),
                                 measurement_std, nside=nside)
true_params = {
    'log_orbital_period':np.log(p_orbit),
    'log_rotation_period':np.log(p_rotation),
    'logit_cos_inc':logit(np.cos(inclination)),
    'logit_cos_obl':logit(np.cos(obliquity)),
    'logit_phi_orb':logit(phi_orb, low=0, high=2*np.pi),
    'logit_obl_orientation':logit(phi_rot, low=0, high=2*np.pi)
}


truth.fix_params(true_params)
p = np.concatenate([np.zeros(truth.nparams), one_point_map])
numeric_lightcurve=truth.lightcurve(p)



# =============ANALYTIC LIGHTCURVE==============================
t0 = t.clock()
analytic_lc = analytic_lightcurve(pt_lat = 0, pt_lng = 0, lat=lat, lng=lng, w_rot=w_rot,
                                          w_orb=w_orb, obq=obliquity, i=inclination,
                                          sol_phase = sol_phase, times=times, nside=nside)
# analytic_lc = analytic_lightcurve(lat=lat, lng=lng, w_rot=w_rot,
#                                           w_orb=w_orb, obq=obliquity, i=inclination,
#                                           sol_phase = sol_phase, times=times, nside=nside)
t1 = t.clock()

# =======================RUN TIMES===============================
analytic_time = t1-t0

t2 = t.clock()

# Numeric solution:
truth = IlluminationMapPosterior(times, np.zeros_like(times),
                                 measurement_std, nside=nside)
truth.fix_params(true_params)
p = np.concatenate([np.zeros(truth.nparams), one_point_map])
numeric =truth.lightcurve(p)
t3 = t.clock()
numeric_time = t3-t2

print("Run times: ")
print("Numeric time: %f. Analytic time: %f" %(numeric_time, analytic_time))


def plot_curves(analytic_lc, numeric_lightcurve, times, savefig):

    # First, obtaining the residuals
    residuals = numeric_lightcurve-analytic_lc

    # Then, plotting the two lightcurves on the same plot
    f, (ax1, ax2) = plt.subplots(2, sharex=True, sharey=False)
    f.suptitle("Analytic vs. Numerical Lightcurves")
    ax1.plot(times, analytic_lc, label='analytic lightcurve')
    ax1.plot(times, numeric_lightcurve, label='numeric lightcurve')
    ax1.legend()
    ax1.set_ylabel('reflectance')
    plt.xlabel('time(h)')

    ax2.plot(times, residuals, label='residuals')
    ax2.set_ylabel('residuals')
    ax2.legend()
    ax2.set_ylim(-0.02, 0.02)
    f.subplots_adjust(hspace=0)
    # Save the fig if true
    if savefig==True:
        plt.savefig('optimal_analytic.pdf')
    plt.show()



# Run times for an entire map:
numeric_time = t3-t2
analytic_time = 0
longitudes = []
latitudes = []
# First, need to convert each entry in the healpix map into lat, lng
for i in range (0, one_point_map.size):
    latitudes.append(map.healpix_to_lonlat([i])[0])
    longitudes.append(map.healpix_to_lonlat([i])[1])


t0=t.clock()
for i in range (0, one_point_map.size):
    analytic_lc_time = analytic_lightcurve(pt_lat=latitudes[i], pt_lng=longitudes[i], lat=lat, lng=lng, w_rot=w_rot,
                                              w_orb=w_orb, obq=obliquity, i=inclination,
                                              sol_phase=sol_phase, times=times,nside=nside)

    # If the map points = the bright spot, then display the lightcurve plot. Else, the lightcurve will be 0 so there
    # will be no plot.
    if i == bright_spot:
        plot_curves(analytic_lc_time, numeric_lightcurve, times, savefig=False)

t1=t.clock()
analytic_time = t1-t0


# PRINTING THE RESULTS FOR THE ENTIRE MAP
print("---------Run-time for a whole map (192 pixels)--------------")
print("Analytic time: %f" %analytic_time)
print("Numeric time: %f" %numeric_time)
print("Residuals: ")
residuals = np.sum(abs(np.asarray(numeric_lightcurve)-np.asarray(analytic_lc)))
print(residuals)




