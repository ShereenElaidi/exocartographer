from analytic_lc import analytic_lightcurve_2
import numpy as np
import healpy as hp
from exocartographer.gp_map import draw_map
from exocartographer import IlluminationMapPosterior
from exocartographer.util import logit, inv_logit
from astropy_healpix import HEALPix
import time as t

# Helper function
def numeric_sln(simulated_map, times, measurement_std, nside, p_orbit, p_rotation, inclination, obliquity, phi_orb, phi_rot):
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
    p = np.concatenate([np.zeros(truth.nparams), simulated_map])
    numeric_lightcurve = truth.lightcurve(p)
    return numeric_lightcurve


def optimize(bright_spot, simulated_map):
    nside = 4
    map = HEALPix(nside=nside, order='nested')

    #Inputs
    p_rotation = 1
    p_orbit = 365.256363
    times = np.linspace(start=0, stop=365, num=1400)
    measurement_std = 0.001
    lat, lng = map.healpix_to_lonlat([bright_spot])
    w_rot = 2*np.pi/p_rotation
    w_orb = 2*np.pi/p_orbit
    inclination = np.pi/2
    obliquity = 0
    sol_phase = 1
    p_rotation = 2*np.pi/w_rot
    p_orb = 2*np.pi/w_orb
    phi_orb = abs(1*(np.pi/2) + sol_phase)       # Unknown parameter 1
    phi_rot = abs(1*(np.pi/2) - sol_phase)       # Unknown parameter 2

    # Analytic lightcurve
    analytic = analytic_lightcurve_2(lat=lat, lng=lng, w_rot=w_rot,
                                              w_orb=w_orb, obq=obliquity, i=inclination,
                                              sol_phase = sol_phase, times=times, nside=nside)
    # Potential values
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

    # Initializing the variables that we need
    opt_residual = 1000
    opt_phi_orb = 0
    opt_phi_rot = 0
    sgn_1 = '+'
    sgn_2 = '-'

    # Optimizng for this single bright spot

    for i in range(0, pot_phi_orb.size):
        for j in range(0, pot_phi_rot.size):
            print("Working...On loop i=%f, j=%f" %(i, j))
            curr_phi_orb = pot_phi_orb.item(i)
            curr_phi_rot = pot_phi_rot.item(j)

            # First: (++)
            phi_orb = abs(curr_phi_orb + sol_phase)
            phi_rot = abs(curr_phi_rot + sol_phase)

            # Numeric lightcurve computation
            numeric = numeric_sln(simulated_map, times, measurement_std, nside, p_orbit, p_rotation,
                                  inclination, obliquity, phi_orb, phi_rot)

            # Resdiuals
            residuals = numeric - analytic
            residuals = np.abs(residuals)
            curr_residuals = np.sum(residuals)

            if curr_residuals < opt_residual:
                opt_residual = curr_residuals
                opt_phi_orb = curr_phi_orb
                opt_phi_rot = curr_phi_rot
                sgn_1 = '+'
                sgn_2 = '+'

            # Second: (--)
            phi_orb = abs(curr_phi_orb - sol_phase)
            phi_rot = abs(curr_phi_rot - sol_phase)

            # Numeric ligtcurve compuation
            numeric = numeric_sln(simulated_map, times, measurement_std, nside, p_orbit, p_rotation,
                                  inclination, obliquity, phi_orb, phi_rot)

            # Resdiuals
            residuals = numeric - analytic
            residuals = np.abs(residuals)
            curr_residuals = np.sum(residuals)

            if curr_residuals < opt_residual:
                opt_residual = curr_residuals
                opt_phi_orb = curr_phi_orb
                opt_phi_rot = curr_phi_rot
                sgn_1 = '-'
                sgn_2 = '-'

            # Third: (-,+)
            phi_orb = abs(curr_phi_orb - sol_phase)
            phi_rot = abs(curr_phi_rot + sol_phase)

            # Numeric ligtcurve compuation
            numeric = numeric_sln(simulated_map, times, measurement_std, nside, p_orbit, p_rotation,
                                  inclination, obliquity, phi_orb, phi_rot)

            # Resdiuals
            residuals = numeric - analytic
            residuals = np.abs(residuals)
            curr_residuals = np.sum(residuals)

            if curr_residuals < opt_residual:
                opt_residual = curr_residuals
                opt_phi_orb = curr_phi_orb
                opt_phi_rot = curr_phi_rot
                sgn_1 = '-'
                sgn_2 = '+'

            # Fourth: (+, -)
            phi_orb = abs(curr_phi_orb + sol_phase)
            phi_rot = abs(curr_phi_rot - sol_phase)

            # Numeric ligtcurve compuation
            numeric = numeric_sln(simulated_map, times, measurement_std, nside, p_orbit, p_rotation,
                                  inclination, obliquity, phi_orb, phi_rot)

            # Resdiuals
            residuals = numeric - analytic
            residuals = np.abs(residuals)
            curr_residuals = np.sum(residuals)

            if curr_residuals < opt_residual:
                opt_residual = curr_residuals
                opt_phi_orb = curr_phi_orb
                opt_phi_rot = curr_phi_rot
                sgn_1 = '+'
                sgn_2 = '-'

    print("Results for i=%f:" % i)
    print("Optimal phi_orb: %f" % opt_phi_orb)
    print("Optimal phi_rot: %f" % opt_phi_rot)
    print("Optimal residuals: %f" % opt_residual)
    print("Sgn 1 %s" % sgn_1)
    print("Sgn 2 %s" % sgn_2)
    return opt_residual



opt_residuals = 1000000000
nside = 4
whitenoise_relative_amp = 0.02
length_scale = 30. * np.pi / 180
albedo_mean = .5
albedo_std = 0.2
while True:
    simulated_map = draw_map(nside, albedo_mean,
                             albedo_std, whitenoise_relative_amp,
                             length_scale)
    if min(simulated_map) > 0 and max(simulated_map) < 1:
        break


opt_k = 0

for k in range (0, simulated_map.size):
    bright_spot = k
    nside = 4
    whitenoise_relative_amp = 0.02
    length_scale = 30. * np.pi / 180
    albedo_mean = .5
    albedo_std = 0.2
    while True:
        simulated_map = draw_map(nside, albedo_mean,
                                 albedo_std, whitenoise_relative_amp,
                                 length_scale)
        if min(simulated_map) > 0 and max(simulated_map) < 1:
            break
    for i in range(0, simulated_map.size):
        if i == bright_spot:
            simulated_map[i] = 1.00
        else:
            simulated_map[i] = 0.00

    residuals = optimize(k, simulated_map)
    print(residuals)


    if residuals < opt_residuals:
        opt_residuals=residuals
        opt_k = k
        print("New optimum found at k = %f" % bright_spot)


print("Results of optimization: ")
print("k=%f" % k)
print("Optimal residuals:")
print(opt_residuals)















