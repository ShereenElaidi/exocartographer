from analytic_lc import analytic_lightcurve_2
import numpy as np
from exocartographer.gp_map import draw_map
from exocartographer import IlluminationMapPosterior
from exocartographer.util import logit, inv_logit
from astropy_healpix import HEALPix
import time as t
from astropy.coordinates import Angle
import astropy.units as u


t0 = t.clock()
# Map
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
bright_spot = 12
for i in range(0, one_point_map.size):
    if i == bright_spot:
        one_point_map[i] = 1.00
    else:
        one_point_map[i] = 0.00
map = HEALPix(nside=nside, order='ring')

#Inputs
p_rotation = 1
p_orbit = 365.256363
times = np.linspace(start=0.0, stop=10.0, num=1400)
measurement_std = 0.001
lat, lng = map.healpix_to_lonlat([bright_spot])
a = Angle((np.pi/2)*u.rad)
colat = a-lat
w_rot = 2*np.pi/p_rotation
w_orb = 2*np.pi/p_orbit
inclination = 0
obliquity = 0
sol_phase = np.pi/6
# p_rotation = 2*np.pi/w_rot
# p_orb = 2*np.pi/w_orb
phi_orb = abs(1*(np.pi/2) + sol_phase)       # Unknown parameter 1
phi_rot = abs(1*(np.pi/2) - sol_phase)       # Unknown parameter 2

# Numeric lightcurve
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
numeric_lightcurve = truth.lightcurve(p)

# Analytic lightcurve
analytic_lightcurve = analytic_lightcurve_2(colat=colat, lng=lng, w_rot=w_rot,
                                          w_orb=w_orb, obq=obliquity, i=inclination,
                                          sol_phase = sol_phase, times=times, nside=nside)

# Optimization part of the script
pot_phi_orb = np.array([0, np.pi/6, np.pi/4, np.pi/3, np.pi/2,
                        2*(np.pi/3), 3*(np.pi/4), 5*(np.pi/6),
                        np.pi, 7*(np.pi/6), 5*(np.pi/4),
                        4*(np.pi/3), 3*(np.pi/2), 5*(np.pi/3),
                        7*(np.pi/4), 11*(np.pi/6), 2*np.pi])
pot_phi_rot = np.array([0, np.pi/6, np.pi/4, np.pi/3, np.pi/2,
                        2*(np.pi/3), 3*(np.pi/4), 5*(np.pi/6),
                        np.pi, 7*(np.pi/6), 5*(np.pi/4),
                        4*(np.pi/3), 3*(np.pi/2), 5*(np.pi/3),
                        7*(np.pi/4), 11*(np.pi/6), 2*np.pi])
opt_residual = 1000
opt_phi_orb = 0
opt_phi_rot = 0
sgn_1 = '+'
sgn_2 = '-'

for i in range(0, pot_phi_orb.size):
    for j in range(0, pot_phi_rot.size):
        print("Working... On loop i = %f j= %f"%(i,j))
        curr_phi_orb = pot_phi_orb.item(i)
        curr_phi_rot = pot_phi_rot.item(j)
        # First: ++
        phi_orb = abs(curr_phi_orb+sol_phase)
        phi_rot = abs(curr_phi_rot+sol_phase)

        # Numeric lightcurve computation:
        truth = IlluminationMapPosterior(times, np.zeros_like(times),
                                         measurement_std, nside=nside)
        true_params = {
            'log_orbital_period': np.log(p_orbit),
            'log_rotation_period': np.log(p_rotation),
            'logit_cos_inc': logit(np.cos(inclination)),
            'logit_cos_obl': logit(np.cos(obliquity)),
            'logit_phi_orb': logit(phi_orb, low=0, high=2 * np.pi),
            'logit_obl_orientation': logit(phi_rot, low=0, high=2 * np.pi)
        }
        truth.fix_params(true_params)
        p = np.concatenate([np.zeros(truth.nparams), one_point_map])
        numeric_lightcurve = truth.lightcurve(p)

        # Residuals
        residuals = numeric_lightcurve-analytic_lightcurve
        residuals = np.absolute(residuals)
        curr_residuals = np.sum(residuals)

        if curr_residuals < opt_residual:
            opt_residual=curr_residuals
            opt_phi_orb = curr_phi_orb
            opt_phi_rot = curr_phi_rot
            sgn_1 = '+'
            sgn_2 = '+'


        # Second: --

        phi_orb = abs(curr_phi_orb-sol_phase)
        phi_rot = abs(curr_phi_rot-sol_phase)
        # Numeric lightcurve computation:
        truth = IlluminationMapPosterior(times, np.zeros_like(times),
                                         measurement_std, nside=nside)
        true_params = {
            'log_orbital_period': np.log(p_orbit),
            'log_rotation_period': np.log(p_rotation),
            'logit_cos_inc': logit(np.cos(inclination)),
            'logit_cos_obl': logit(np.cos(obliquity)),
            'logit_phi_orb': logit(phi_orb, low=0, high=2 * np.pi),
            'logit_obl_orientation': logit(phi_rot, low=0, high=2 * np.pi)
        }
        truth.fix_params(true_params)
        p = np.concatenate([np.zeros(truth.nparams), one_point_map])
        numeric_lightcurve = truth.lightcurve(p)

        # Residuals
        residuals = numeric_lightcurve-analytic_lightcurve
        residuals = np.absolute(residuals)
        curr_residuals = np.sum(residuals)

        if curr_residuals < opt_residual:
            opt_residual=curr_residuals
            opt_phi_orb = curr_phi_orb
            opt_phi_rot = curr_phi_rot
            sgn_1 = '-'
            sgn_2 = '-'

        # Third: +-
        phi_orb = abs(curr_phi_orb+sol_phase)
        phi_rot = abs(curr_phi_rot-sol_phase)
        # Numeric lightcurve computation:
        truth = IlluminationMapPosterior(times, np.zeros_like(times),
                                         measurement_std, nside=nside)
        true_params = {
            'log_orbital_period': np.log(p_orbit),
            'log_rotation_period': np.log(p_rotation),
            'logit_cos_inc': logit(np.cos(inclination)),
            'logit_cos_obl': logit(np.cos(obliquity)),
            'logit_phi_orb': logit(phi_orb, low=0, high=2 * np.pi),
            'logit_obl_orientation': logit(phi_rot, low=0, high=2 * np.pi)
        }
        truth.fix_params(true_params)
        p = np.concatenate([np.zeros(truth.nparams), one_point_map])
        numeric_lightcurve = truth.lightcurve(p)

        # Residuals
        residuals = numeric_lightcurve-analytic_lightcurve
        residuals = np.absolute(residuals)
        curr_residuals = np.sum(residuals)

        if curr_residuals < opt_residual:
            opt_residual=curr_residuals
            opt_phi_orb = curr_phi_orb
            opt_phi_rot = curr_phi_rot
            sgn_1 = '+'
            sgn_2 = '-'

        # Fourth: -+
        phi_orb = abs(curr_phi_orb-sol_phase)
        phi_rot = abs(curr_phi_rot+sol_phase)
        # Numeric lightcurve computation:
        truth = IlluminationMapPosterior(times, np.zeros_like(times),
                                         measurement_std, nside=nside)
        true_params = {
            'log_orbital_period': np.log(p_orbit),
            'log_rotation_period': np.log(p_rotation),
            'logit_cos_inc': logit(np.cos(inclination)),
            'logit_cos_obl': logit(np.cos(obliquity)),
            'logit_phi_orb': logit(phi_orb, low=0, high=2 * np.pi),
            'logit_obl_orientation': logit(phi_rot, low=0, high=2 * np.pi)
        }
        truth.fix_params(true_params)
        p = np.concatenate([np.zeros(truth.nparams), one_point_map])
        numeric_lightcurve = truth.lightcurve(p)

        # Residuals
        residuals = numeric_lightcurve-analytic_lightcurve
        residuals = np.absolute(residuals)
        curr_residuals = np.sum(residuals)

        if curr_residuals < opt_residual:
            opt_residual=curr_residuals
            opt_phi_orb = curr_phi_orb
            opt_phi_rot = curr_phi_rot
            sgn_1 = '-'
            sgn_2 = '+'
t1= t.clock()
run_time = t1-t0
print("Running time: %f" % run_time)
print("Results:")
print("Optimal phi_orb: %f" % opt_phi_orb )
print("Optimal phi_rot: %f" % opt_phi_rot)
print("Optimal residuals: %f" % opt_residual)
print("Sgn 1: %s" % sgn_1)
print("Sgn 2: %s" % sgn_2)