import matplotlib.pyplot as plt
from analytic_lc import analytic_lightcurve
import numpy as np
import healpy as hp
from scipy.optimize import minimize
from exocartographer.gp_map import draw_map
from exocartographer import IlluminationMapPosterior
from exocartographer.util import logit, inv_logit
import pdb
from astropy_healpix import HEALPix
# resolution
nside = 4

# Gaussian process properties
whitenoise_relative_amp = 0.02
length_scale = 30.*np.pi/2
albedo_mean = 0.5
albedo_std = 0.2

# drawing an albedo map
while True:
    simulated_map = draw_map(nside, albedo_mean, albedo_std,
                             whitenoise_relative_amp, length_scale)
    if min(simulated_map) > 0 and max(simulated_map) < 1:
        break
#
# # simulating a map that's all black everywhere except for one bright spot
#
for i in range(0, simulated_map.size):
    if i == 100:
        simulated_map[i] = 1.00
    else:
        simulated_map[i] = 0.000
hp.mollview(simulated_map, title='albedo', cmap='gist_gray')

# obtaining the lon and lat of this point
map  = HEALPix(nside=nside, order='nested')
pt_lng, pt_lat = map.healpix_to_lonlat([100])
print("lng")
print(pt_lng)
print("lat")
print(pt_lat)

# Set orbital properties
p_rotation = 23.934
p_orbit = 365.256363 * 24.0
phi_orb = np.pi
# inclination = np.pi/2
# obliquity = 90. * np.pi/180.0
obliquity = np.pi/4
inclination = np.pi/2
phi_rot = np.pi

# Observation schedule
cadence = p_rotation/4.
epoch_duration = p_orbit
times = np.linspace(0, epoch_duration, epoch_duration/cadence)

times = np.arange(1.0, 20.0, 0.1)
times = times/p_rotation

# Measurement uncertainties
measurement_std = 0.001

# Use a posterior instance for easy lightcurve generation
truth = IlluminationMapPosterior(times, np.zeros_like(times),
                                 measurement_std, nside=nside)

true_params = {
    'log_orbital_period':np.log(p_orbit),
    'log_rotation_period':np.log(p_rotation),
    'logit_cos_inc':logit(np.cos(inclination)),
    'logit_cos_obl':logit(np.cos(obliquity)),
    'logit_phi_orb':logit(phi_orb, low=0, high=2*np.pi),
    'logit_obl_orientation':logit(phi_rot, low=0, high=2*np.pi)}
truth.fix_params(true_params)

p = np.concatenate([np.zeros(truth.nparams), simulated_map])
lightcurve = truth.lightcurve(p)

plt.figure(figsize=(16,6))
plt.subplot(2, 1, 1)
# plt.figure(figsize=(16, 3))
plt.plot(times, lightcurve, lw=0.5)
plt.xlim(0, p_orbit/p_rotation)
plt.ylim(ymin=0)
plt.xlabel(r'time$/P_\mathrm{rot}$')
plt.ylabel('reflectance')



w_rot = 1/p_rotation
w_orb = 1/p_orbit
time_data = np.arange(0.0, 20.0, 0.1)
# time_data = times
obq = 0
i = np.pi/4
sol_phase = 1
lng = 3.53429
lat = np.pi/2--0.167448

flux = analytic_lightcurve(lat, lng, w_rot, w_orb, obq, i, sol_phase, time_data)
print(flux)

plt.subplot(2, 1, 2)
# plt.figure(figsize=(16, 3))
plt.title("Analytic lightcurve")
plt.plot(time_data,flux)
plt.ylabel('flux')
plt.xlabel(r'time$/P_\mathrm{rot}$')
plt.show()

