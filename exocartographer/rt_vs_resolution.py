# Import statements
import time as t
import numpy as np
import matplotlib.pyplot as plt
from analytic_lc import analytic_lightcurve
import healpy as hp
from exocartographer.gp_map import draw_map
from exocartographer import IlluminationMapPosterior
from exocartographer.util import logit, inv_logit
from astropy_healpix import HEALPix
import pdb


# Potential Nside values = 1,2,4,8

# Map:
nside = 1
map = HEALPix(nside=nside, order='nested')
whitenoise_relative_amp = 0.2
length_scale = 30*np.pi/2
albedo_mean = 0.5
albedo_std = 0.2
while True:
    one_point_map = draw_map(nside, albedo_mean, albedo_std,
                             whitenoise_relative_amp, length_scale)
    if min(one_point_map) > 0 and max(one_point_map) <1:
        break
bright_spot = 1
for i in range(0, one_point_map.size):
    if i == bright_spot:
        one_point_map[i] = 1.00
    else:
        one_point_map[i] = 0.00

# Parameters
p_rotation = 1
p_orbit = 365.256363
lat, lng = map.healpix_to_lonlat([bright_spot])
w_rot = 2*np.pi/p_rotation
w_orb = 2*np.pi/p_orbit
inclination = np.pi/2
obliquity = 0
sol_phase = 1
p_rotation = 2*np.pi/w_rot
p_orb = 2*np.pi/w_orb
times = np.linspace(start=0.0, stop=1.0, num=1400)
measurement_std = 0.001

# DO NOT CHANGE THIS.
phi_orb = abs(5*(np.pi/3) + sol_phase)
phi_rot = abs(2*np.pi - sol_phase)


# NUMERIC
true_params = {
    'log_orbital_period':np.log(p_orbit),
    'log_rotation_period':np.log(p_rotation),
    'logit_cos_inc':logit(np.cos(inclination)),
    'logit_cos_obl':logit(np.cos(obliquity)),
    'logit_phi_orb':logit(phi_orb, low=0, high=2*np.pi),
    'logit_obl_orientation':logit(phi_rot, low=0, high=2*np.pi)
}
truth = IlluminationMapPosterior(times, np.zeros_like(times),
                                 measurement_std, nside=nside)
truth.fix_params(true_params)
p = np.concatenate([np.zeros(truth.nparams), one_point_map])
numeric_lightcurve=truth.lightcurve(p)


run_times_analytic = np.array([])
run_times_numeric = np.array([])
nside_resolutions  = np.array([1,2,4,8])
print("Test:")
for i in np.nditer(nside_resolutions):
    nside = i
    whitenoise_relative_amp = 0.2
    length_scale = 30 * np.pi / 2
    albedo_mean = 0.5
    albedo_std = 0.2

    while True:
        one_point_map = draw_map(nside, albedo_mean, albedo_std,
                                 whitenoise_relative_amp, length_scale)
        if min(one_point_map) > 0 and max(one_point_map) < 1:
            break
    bright_spot = 1
    for i in range(0, one_point_map.size):
        if i == bright_spot:
            one_point_map[i] = 1.00
        else:
            one_point_map[i] = 0.00
    map = HEALPix(nside=nside, order='nested')

    # Put numeric code here:
    pdb.set_trace()
    t0 = t.clock()

    truth = IlluminationMapPosterior(times, np.zeros_like(times),
                                     measurement_std, nside=nside)
    truth.fix_params(true_params)
    p = np.concatenate([np.zeros(truth.nparams), one_point_map])
    numeric_lightcurve = truth.lightcurve(p)
    t1 = t.clock()

    num_time = t1-t0
    np.append(run_times_numeric, num_time)

    # Put analytic code here:
    analytic_time = 0
    longitudes = []
    latitudes = []
    for i in range(0, one_point_map.size):
        latitudes.append(map.healpix_to_lonlat([i])[0])
        longitudes.append(map.healpix_to_lonlat([i])[1])

    t2 = t.clock()

    for i in range(0, one_point_map.size):
        analytic_lc_time = analytic_lightcurve(pt_lat=latitudes[i], pt_lng=longitudes[i], lat=lat, lng=lng, w_rot=w_rot,
                                               w_orb=w_orb, obq=obliquity, i=inclination,
                                               sol_phase=sol_phase, times=times)

    t3 = t.clock()
    an_time = t3-t2
    np.append(run_times_analytic, an_time)

print("Results:")
print("Analytic running times:")
print(run_times_analytic)
print("Numeric running times:")
print(run_times_numeric)