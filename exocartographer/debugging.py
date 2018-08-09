# Import statements
import numpy as np
import matplotlib.pyplot as plt
from analytic_lc import analytic_lightcurve_2
import healpy as hp
from exocartographer.gp_map import draw_map
from exocartographer import IlluminationMapPosterior
from exocartographer.util import logit, inv_logit
from astropy_healpix import HEALPix
from astropy.coordinates import Angle
import astropy.units as u
import matplotlib.gridspec  as gridspec
from exocartographer import viewgeom, kernel
import math
from more_debugging import get_visibility_illum, get_points

# Helper Functions

# (1) A function to plot curves
def plot_curves(analytic_lc, numeric_lightcurve, times, savefig):

    # First, obtaining the residuals
    residuals = numeric_lightcurve-analytic_lc

    # Then, plotting the two lightcurves on the same plot
    f, (ax1, ax2) = plt.subplots(2, sharex=True, sharey=True)
    ax1.plot(times, analytic_lc, label='analytic lightcurve')
    ax1.plot(times, numeric_lightcurve, label='numeric lightcurve')
    ax1.legend()
    plt.ylabel('reflectance')
    plt.xlabel('time(days)')

    ax2.plot(times, residuals, label='residuals')
    ax2.legend()
    f.subplots_adjust(hspace=0)
    # Save the fig if true
    if savefig==True:
        plt.savefig('optimal_analytic.pdf')
    plt.show()

# (2) A function to draw the HEALPix map
def draw_albedo_map(nside, view_map, bright_spot):
    nside = nside # map resolution

    # Gaussian process properties
    whitenoise_relative_amp = 0.02
    length_scale = 30. * np.pi / 180
    albedo_mean = .5
    albedo_std = 0.2

    # Draw a valid albedo map (i.e., 0 < albedo < 1)
    while True:
        simulated_map = draw_map(nside, albedo_mean, albedo_std,
                                 whitenoise_relative_amp, length_scale)
        if min(simulated_map) > 0 and max(simulated_map) < 1:
            break

    for i in range(0, simulated_map.size):
        if i == bright_spot:
            simulated_map[i] = 1.00
        else:
            simulated_map[i] = 0.00
    if view_map==True:
        hp.mollview(simulated_map, title='albedo', cmap='gist_gray')
    return simulated_map

# (3) A function to produce the *visibility* plot. This is lifted from the *analytic* code.
def Vnz (colat, lng, i, obq, sol_phase, times):
    times = times
    p_rotation = 23.934
    p_orbit = 365.256363 * 24.0
    w_rot = 2 * np.pi / p_rotation
    w_orb = 2 * np.pi / p_orbit
    i = i
    obq = obq
    sol_phase = sol_phase

    # ------------------
    theta = colat
    phi = lng

    # Appendix C: Derivation of the Time-Varying Geometry
    # (C1)
    cos_theta_knot = (np.cos(i) * np.cos(obq)) + (np.sin(i) * np.sin(obq) * np.cos(sol_phase))
    sin_theta_knot = np.sqrt(1 - (cos_theta_knot ** 2))

    # (C2)
    cos_phi_knot = np.cos(w_rot * times)
    sin_phi_knot = -np.sin(w_rot * times)

    # (C3)
    cos_theta_s = np.sin(obq) * np.cos((w_orb * times) - sol_phase)
    sin_theta_s = np.sqrt(1 - (((np.sin(obq)) ** 2) * ((np.cos((w_orb * times) - sol_phase)) ** 2)))
    cos_phi_s = 0
    sin_phi_s = 0

    if i == np.pi / 2:
        # (C1) Orbital Inclination
        cos_theta_knot = np.sin(obq) * np.cos(sol_phase)

        # (C12)
        cos_phi_s_1 = (np.cos(w_rot * times)) * (
                np.cos(w_orb * times) - cos_theta_knot * np.sin(obq) * np.cos(w_orb * times - sol_phase))
        cos_phi_s_2 = np.sin(w_rot * times) * np.sin(w_orb * times) * np.cos(obq)
        cos_phi_s_3 = sin_theta_knot * (
            np.sqrt(1 - ((((np.sin(obq)) ** 2) * ((np.cos(w_orb * times - sol_phase)) ** 2)))))
        cos_phi_s = (cos_phi_s_1 + cos_phi_s_2) / (cos_phi_s_3)

        # (C13)
        sin_phi_s_1 = -np.sin(w_rot * times) * (
                np.cos(w_orb * times) - cos_theta_knot * np.sin(obq) * np.cos(w_orb * times - sol_phase))
        sin_phi_s_2 = np.cos(w_rot * times) * np.sin(w_orb * times) * np.cos(obq)
        sin_phi_s_3 = sin_theta_knot * (
            np.sqrt(1 - (((np.sin(obq)) ** 2) * ((np.cos(w_orb * times - sol_phase)) ** 2))))
        sin_phi_s = (sin_phi_s_1 + sin_phi_s_2) / sin_phi_s_3

    elif i == 0:
        cos_theta_knot = np.cos(obq)

        # (C14)
        cos_phi_s_1 = -np.cos(w_rot * times) * np.cos(obq) * np.cos(w_orb * times - sol_phase)
        cos_phi_s_2 = np.sin(w_rot * times) * np.sin(w_orb * times - sol_phase)
        cos_phi_s_3 = np.sqrt(1 - ((((np.sin(obq)) ** 2)) * ((np.cos(w_orb * times - sol_phase)) ** 2)))
        cos_phi_s = (cos_phi_s_1 - cos_phi_s_2) / cos_phi_s_3

        # (C15)
        sin_phi_s_1 = np.sin(w_rot * times) * np.cos(obq) * np.cos(w_orb * times - sol_phase)
        sin_phi_s_2 = np.cos(w_rot * times) * np.sin(w_orb * times - sol_phase)
        sin_phi_s_3 = np.sqrt(1 - (((np.sin(obq)) ** 2) * ((np.cos(w_orb * times - sol_phase)) ** 2)))
        sin_phi_s = (sin_phi_s_1 - sin_phi_s_2) / sin_phi_s_3

    elif obq == 0:

        # (C3) Planetary Obliquity
        theta_knot = i
        theta_s = np.pi / 2
        phi_s = (w_orb - w_rot) * times
        cos_theta_knot = np.cos(i)
        cos_theta_s = 0
        cos_phi_s = np.cos((w_orb - w_rot) * times)
        sin_phi_s = np.sin((w_orb - w_rot) * times)

    elif obq == np.pi / 2:
        # (C20)
        cos_theta_knot = np.sin(i) * np.sin(obq) * np.cos(sol_phase)
        sin_theta_knot = np.sqrt(1 - (cos_theta_knot ** 2))

        # (C21)
        cos_theta_s = np.cos(w_orb * times - sol_phase)
        sin_theta_s = np.sin(w_orb * times - sol_phase)

        # (C22)
        cos_phi_s_1 = np.cos(w_rot * times) * (
                np.sin(i) * np.cos(w_orb * times) - cos_theta_knot * np.cos(w_orb * times - sol_phase))
        cos_phi_s_2 = np.sin(w_rot * times) * np.cos(i) * np.sin(w_orb * times - sol_phase)
        cos_phi_s_3 = np.sqrt(1 - (((np.sin(i)) ** 2) * ((np.cos(sol_phase)) ** 2)))
        cos_phi_s_4 = np.sin(w_orb * times - sol_phase)
        cos_phi_s = (cos_phi_s_1 - cos_phi_s_2) / (cos_phi_s_3 * cos_phi_s_4)

        # (C23)
        sin_phi_s_1 = -np.sin(w_rot * times) * (
                np.sin(i) * np.cos(w_orb * times) - cos_theta_knot * np.cos(w_orb * times - sol_phase))
        sin_phi_s_2 = np.cos(w_rot * times) * np.cos(i) * np.sin(w_orb * times - sol_phase)
        sin_phi_s_3 = np.sqrt(1 - (((np.sin(i)) ** 2) * ((np.cos(sol_phase)) ** 2)))
        sin_phi_s_4 = np.sin(w_orb * times - sol_phase)
        sin_phi_s = (sin_phi_s_1 - sin_phi_s_2) / (sin_phi_s_3 * sin_phi_s_4)

    else:
        # (C4)
        cos_phi_s_1 = np.cos(w_rot * times) * (
                np.sin(i) * np.cos(w_orb * times) - cos_theta_knot * np.sin(obq) * np.cos(w_orb * times - sol_phase))
        cos_phi_s_2 = np.sin(w_rot * times) * (
                np.sin(i) * np.sin(w_orb * times) * np.cos(obq) - np.cos(i) * np.sin(obq) * np.sin(
            w_orb * times - sol_phase))
        cos_phi_s_3 = np.sqrt(1 - (np.cos(i) * np.cos(obq) + np.sin(i) * np.sin(obq) * np.cos(sol_phase)) ** 2)
        cos_phi_s_4 = np.sqrt(1 - (((np.sin(obq)) ** 2) * ((np.cos(w_orb * times - sol_phase)) ** 2)))
        cos_phi_s = (cos_phi_s_1 + cos_phi_s_2) / (cos_phi_s_3 * cos_phi_s_4)

        # (C5)
        sin_phi_s_1 = -np.sin(w_rot * times) * (
                np.sin(i) * np.cos(w_orb * times) - cos_theta_knot * np.sin(obq) * np.cos(w_orb * times - sol_phase))
        sin_phi_s_2 = np.cos(w_rot * times) * (
                np.sin(i) * np.sin(w_orb * times) * np.cos(obq) - np.cos(i) * np.sin(obq) * np.sin(
            w_orb * times - sol_phase))
        sin_phi_s_3 = np.sqrt(1 - (np.cos(i) * np.cos(obq) + np.sin(i) * np.sin(obq) * np.cos(sol_phase)) ** 2)
        sin_phi_s_4 = np.sqrt(1 - (((np.sin(obq)) ** 2) * ((np.cos(w_orb * times - sol_phase)) ** 2)))
        sin_phi_s = (sin_phi_s_1 + sin_phi_s_2) / (sin_phi_s_3 * sin_phi_s_4)

    # Delta Maps
    cos_phi_phi_s = np.cos(phi) * cos_phi_s + np.sin(phi) * sin_phi_s
    cos_phi_phi_knot = np.cos(phi) * cos_phi_knot + np.sin(phi) * sin_phi_knot

    # (4)
    Inz = np.sin(theta) * sin_theta_s * cos_phi_phi_s + np.cos(theta) * cos_theta_s

    # (5)
    Vnz = np.sin(theta) * sin_theta_knot * cos_phi_phi_knot + np.cos(theta) * cos_theta_knot
    return np.maximum(Vnz, 0)

# (4) A function to produce the *illumination* plot. This is lifted from the *analytic* code.
def Inz (colat, lng, i, obq, sol_phase, times):
    times = times
    p_rotation = 23.934
    p_orbit = 365.256363 * 24.0
    w_rot = 2 * np.pi / p_rotation
    w_orb = 2 * np.pi / p_orbit
    i = i
    obq = obq
    sol_phase = sol_phase

    # ------------------
    theta = colat
    phi = lng

    # Appendix C: Derivation of the Time-Varying Geometry
    # (C1)
    cos_theta_knot = (np.cos(i) * np.cos(obq)) + (np.sin(i) * np.sin(obq) * np.cos(sol_phase))
    sin_theta_knot = np.sqrt(1 - (cos_theta_knot ** 2))

    # (C2)
    cos_phi_knot = np.cos(w_rot * times)
    sin_phi_knot = -np.sin(w_rot * times)

    # (C3)
    cos_theta_s = np.sin(obq) * np.cos((w_orb * times) - sol_phase)
    sin_theta_s = np.sqrt(1 - (((np.sin(obq)) ** 2) * ((np.cos((w_orb * times) - sol_phase)) ** 2)))
    cos_phi_s = 0
    sin_phi_s = 0

    if i == np.pi / 2:
        # (C1) Orbital Inclination
        cos_theta_knot = np.sin(obq) * np.cos(sol_phase)

        # (C12)
        cos_phi_s_1 = (np.cos(w_rot * times)) * (
                np.cos(w_orb * times) - cos_theta_knot * np.sin(obq) * np.cos(w_orb * times - sol_phase))
        cos_phi_s_2 = np.sin(w_rot * times) * np.sin(w_orb * times) * np.cos(obq)
        cos_phi_s_3 = sin_theta_knot * (
            np.sqrt(1 - ((((np.sin(obq)) ** 2) * ((np.cos(w_orb * times - sol_phase)) ** 2)))))
        cos_phi_s = (cos_phi_s_1 + cos_phi_s_2) / (cos_phi_s_3)

        # (C13)
        sin_phi_s_1 = -np.sin(w_rot * times) * (
                np.cos(w_orb * times) - cos_theta_knot * np.sin(obq) * np.cos(w_orb * times - sol_phase))
        sin_phi_s_2 = np.cos(w_rot * times) * np.sin(w_orb * times) * np.cos(obq)
        sin_phi_s_3 = sin_theta_knot * (
            np.sqrt(1 - (((np.sin(obq)) ** 2) * ((np.cos(w_orb * times - sol_phase)) ** 2))))
        sin_phi_s = (sin_phi_s_1 + sin_phi_s_2) / sin_phi_s_3

    elif i == 0:
        cos_theta_knot = np.cos(obq)

        # (C14)
        cos_phi_s_1 = -np.cos(w_rot * times) * np.cos(obq) * np.cos(w_orb * times - sol_phase)
        cos_phi_s_2 = np.sin(w_rot * times) * np.sin(w_orb * times - sol_phase)
        cos_phi_s_3 = np.sqrt(1 - ((((np.sin(obq)) ** 2)) * ((np.cos(w_orb * times - sol_phase)) ** 2)))
        cos_phi_s = (cos_phi_s_1 - cos_phi_s_2) / cos_phi_s_3

        # (C15)
        sin_phi_s_1 = np.sin(w_rot * times) * np.cos(obq) * np.cos(w_orb * times - sol_phase)
        sin_phi_s_2 = np.cos(w_rot * times) * np.sin(w_orb * times - sol_phase)
        sin_phi_s_3 = np.sqrt(1 - (((np.sin(obq)) ** 2) * ((np.cos(w_orb * times - sol_phase)) ** 2)))
        sin_phi_s = (sin_phi_s_1 - sin_phi_s_2) / sin_phi_s_3

    elif obq == 0:

        # (C3) Planetary Obliquity
        theta_knot = i
        theta_s = np.pi / 2
        phi_s = (w_orb - w_rot) * times
        cos_theta_knot = np.cos(i)
        cos_theta_s = 0
        cos_phi_s = np.cos((w_orb - w_rot) * times)
        sin_phi_s = np.sin((w_orb - w_rot) * times)

    elif obq == np.pi / 2:
        # (C20)
        cos_theta_knot = np.sin(i) * np.sin(obq) * np.cos(sol_phase)
        sin_theta_knot = np.sqrt(1 - (cos_theta_knot ** 2))

        # (C21)
        cos_theta_s = np.cos(w_orb * times - sol_phase)
        sin_theta_s = np.sin(w_orb * times - sol_phase)

        # (C22)
        cos_phi_s_1 = np.cos(w_rot * times) * (
                np.sin(i) * np.cos(w_orb * times) - cos_theta_knot * np.cos(w_orb * times - sol_phase))
        cos_phi_s_2 = np.sin(w_rot * times) * np.cos(i) * np.sin(w_orb * times - sol_phase)
        cos_phi_s_3 = np.sqrt(1 - (((np.sin(i)) ** 2) * ((np.cos(sol_phase)) ** 2)))
        cos_phi_s_4 = np.sin(w_orb * times - sol_phase)
        cos_phi_s = (cos_phi_s_1 - cos_phi_s_2) / (cos_phi_s_3 * cos_phi_s_4)

        # (C23)
        sin_phi_s_1 = -np.sin(w_rot * times) * (
                np.sin(i) * np.cos(w_orb * times) - cos_theta_knot * np.cos(w_orb * times - sol_phase))
        sin_phi_s_2 = np.cos(w_rot * times) * np.cos(i) * np.sin(w_orb * times - sol_phase)
        sin_phi_s_3 = np.sqrt(1 - (((np.sin(i)) ** 2) * ((np.cos(sol_phase)) ** 2)))
        sin_phi_s_4 = np.sin(w_orb * times - sol_phase)
        sin_phi_s = (sin_phi_s_1 - sin_phi_s_2) / (sin_phi_s_3 * sin_phi_s_4)

    else:
        # (C4)
        cos_phi_s_1 = np.cos(w_rot * times) * (
                np.sin(i) * np.cos(w_orb * times) - cos_theta_knot * np.sin(obq) * np.cos(w_orb * times - sol_phase))
        cos_phi_s_2 = np.sin(w_rot * times) * (
                np.sin(i) * np.sin(w_orb * times) * np.cos(obq) - np.cos(i) * np.sin(obq) * np.sin(
            w_orb * times - sol_phase))
        cos_phi_s_3 = np.sqrt(1 - (np.cos(i) * np.cos(obq) + np.sin(i) * np.sin(obq) * np.cos(sol_phase)) ** 2)
        cos_phi_s_4 = np.sqrt(1 - (((np.sin(obq)) ** 2) * ((np.cos(w_orb * times - sol_phase)) ** 2)))
        cos_phi_s = (cos_phi_s_1 + cos_phi_s_2) / (cos_phi_s_3 * cos_phi_s_4)

        # (C5)
        sin_phi_s_1 = -np.sin(w_rot * times) * (
                np.sin(i) * np.cos(w_orb * times) - cos_theta_knot * np.sin(obq) * np.cos(w_orb * times - sol_phase))
        sin_phi_s_2 = np.cos(w_rot * times) * (
                np.sin(i) * np.sin(w_orb * times) * np.cos(obq) - np.cos(i) * np.sin(obq) * np.sin(
            w_orb * times - sol_phase))
        sin_phi_s_3 = np.sqrt(1 - (np.cos(i) * np.cos(obq) + np.sin(i) * np.sin(obq) * np.cos(sol_phase)) ** 2)
        sin_phi_s_4 = np.sqrt(1 - (((np.sin(obq)) ** 2) * ((np.cos(w_orb * times - sol_phase)) ** 2)))
        sin_phi_s = (sin_phi_s_1 + sin_phi_s_2) / (sin_phi_s_3 * sin_phi_s_4)

    # Delta Maps
    cos_phi_phi_s = np.cos(phi) * cos_phi_s + np.sin(phi) * sin_phi_s
    cos_phi_phi_knot = np.cos(phi) * cos_phi_knot + np.sin(phi) * sin_phi_knot

    # (4)
    Inz = np.sin(theta) * sin_theta_s * cos_phi_phi_s + np.cos(theta) * cos_theta_s
    return np.maximum(0, Inz)

# Setting the orbital properties for the numeric solution
nside = 4

# From the numeric solution: exocartographer uses order = ring.
order = 'ring'
map = HEALPix(nside=nside, order=order)

# The unit here will be hours
p_rotation = 23.934
p_orbit = 365.256363 * 24.0

# Setting the inputs up for the analytic and numeric solutions

# Creating an evenly-spaced time array containing num # of points
times = np.linspace(start=0.0, stop=24.0, num=1400)
measurement_std = 0.001

# Setting the input parameters
w_rot = 2*np.pi/p_rotation
w_orb = 2*np.pi/p_orbit
inclination = np.pi/4
obliquity = np.pi/6
sol_phase = np.pi/6

# Using the relation from the code here
phi_orb = 0                     # xi_0
obl_orientation = sol_phase     # xi_s

# Drawing the simulated map. Same longitude, but varying co-latitudes
lat_1 = Angle((np.pi/2)*u.rad)
lat_2 = Angle((np.pi/3)*u.rad)
lat_3 = Angle((np.pi/4)*u.rad)
lat_4 = Angle((np.pi/6)*u.rad)
lat_5 = Angle((0)*u.rad)

# Converting the latitude values nto co-latitude values
a = Angle((np.pi/2)*u.rad)
colat_1 = a-lat_1
colat_2 = a-lat_2
colat_3 = a-lat_3
colat_4 = a-lat_4
colat_5 = a-lat_5

# Setting longitude to be -pi/2
one_lng = Angle(-(0.8)*u.rad)

# Obtaining the bright spot for each point
bs_1 = map.lonlat_to_healpix([one_lng]*u.rad, [lat_1]*u.rad)
bs_2 = map.lonlat_to_healpix([one_lng]*u.rad, [lat_2]*u.rad)
bs_3 = map.lonlat_to_healpix([one_lng]*u.rad, [lat_3]*u.rad)
bs_4 = map.lonlat_to_healpix([one_lng]*u.rad, [lat_4]*u.rad)
bs_5 = map.lonlat_to_healpix([one_lng]*u.rad, [lat_5]*u.rad)

# Drawing the simulated map for each bright spot
sim_map_1 = draw_albedo_map(nside=nside, view_map=False, bright_spot=bs_1)
sim_map_2 = draw_albedo_map(nside=nside, view_map=False, bright_spot=bs_2)
sim_map_3 = draw_albedo_map(nside=nside, view_map=False, bright_spot=bs_3)
sim_map_4 = draw_albedo_map(nside=nside, view_map=False, bright_spot=bs_4)
sim_map_5 = draw_albedo_map(nside=nside, view_map=False, bright_spot=bs_5)


# Analytic solution
analytic_1 = analytic_lightcurve_2(colat=colat_1, lng=one_lng, w_rot=w_rot, w_orb = w_orb,
                                   obq = obliquity, i=inclination, sol_phase=sol_phase,
                                   times=times, nside=nside)
analytic_2 = analytic_lightcurve_2(colat=colat_2, lng=one_lng, w_rot=w_rot, w_orb = w_orb,
                                   obq=obliquity, i=inclination, sol_phase=sol_phase,
                                   times=times, nside=nside)
analytic_3 = analytic_lightcurve_2(colat=colat_3, lng=one_lng, w_rot =w_rot, w_orb=w_orb,
                                   obq=obliquity, i=inclination, sol_phase=sol_phase,
                                   times=times, nside=nside)
analytic_4 = analytic_lightcurve_2(colat=colat_4, lng=one_lng, w_rot=w_rot, w_orb=w_orb,
                                   obq=obliquity, i=inclination, sol_phase=sol_phase,
                                   times=times, nside=nside)
analytic_5 = analytic_lightcurve_2(colat=colat_5, lng=one_lng, w_rot=w_rot, w_orb=w_orb,
                                   obq=obliquity, i=inclination, sol_phase=sol_phase,
                                   times=times, nside=nside)

# Numeric solution
truth = IlluminationMapPosterior(times, np.zeros_like(times),
                                 measurement_std, nside=nside)
true_params = {
    'log_orbital_period':np.log(p_orbit),
    'log_rotation_period':np.log(p_rotation),
    'logit_cos_inc':logit(np.cos(inclination)),
    'logit_cos_obl':logit(np.cos(obliquity)),
    'logit_phi_orb':logit(phi_orb, low=0, high=2*np.pi),
    'logit_obl_orientation':logit(obl_orientation, low=0, high=2*np.pi)
}


truth.fix_params(true_params)
p_1 = np.concatenate([np.zeros(truth.nparams), sim_map_1])
p_2 = np.concatenate([np.zeros(truth.nparams), sim_map_2])
p_3 = np.concatenate([np.zeros(truth.nparams), sim_map_3])
p_4 = np.concatenate([np.zeros(truth.nparams), sim_map_4])
p_5 = np.concatenate([np.zeros(truth.nparams), sim_map_5])

numeric_1 = truth.lightcurve(p_1)
numeric_2 = truth.lightcurve(p_2)
numeric_3 = truth.lightcurve(p_3)
numeric_4 = truth.lightcurve(p_4)
numeric_5 = truth.lightcurve(p_5)


# Trying to use the analytic_visibility_illumination kernel
num_1 = truth.analytic_visibility_illumination_matrix(p_1)
num_2 = truth.analytic_visibility_illumination_matrix(p_2)
num_3 = truth.analytic_visibility_illumination_matrix(p_3)
num_4 = truth.analytic_visibility_illumination_matrix(p_4)
num_5 = truth.analytic_visibility_illumination_matrix(p_5)


# Making the big plot
fig=plt.figure()
ax= plt.subplot2grid((5, 2), (0,0))

# Subplot 1
ax1 = plt.subplot2grid((5, 2), (0,0), colspan=1, rowspan=1)
ax1.plot(times, analytic_1, label='analytic 1', color='#f20e0e')
ax1.plot(times, numeric_1, label='numeric 1', color='#000000')
ax1.grid(color='#e6e8ed', linestyle='-', linewidth=2)
ax1.legend()


# Subplot 2
ax2 = plt.subplot2grid((5, 2), (1,0), colspan=1, rowspan=1, sharex=ax1, sharey=ax1)
ax2.plot(times, analytic_2, label='analytic 2', color='#ff7f00')
ax2.plot(times, numeric_2, label='numeric 2', color='#000000')
ax2.grid(color='#e6e8ed', linestyle='-', linewidth=2)
ax2.legend()

# Subplot 3
ax3 = plt.subplot2grid((5, 2), (2,0), colspan=1, rowspan=1, sharex=ax1, sharey=ax1)
ax3.plot(times, analytic_3, label='analytic 3', color='#ffc300')
ax3.plot(times, numeric_3, label='numeric 3', color='#000000')
ax3.grid(color='#e6e8ed', linestyle='-', linewidth=2)
ax3.legend()
ax3.set_ylabel('reflectance')

# Subplot 4
ax4 = plt.subplot2grid((5, 2), (3,0), colspan=1, rowspan=1, sharex=ax1, sharey=ax1)
ax4.plot(times, analytic_4, label='analytic 4', color='#03bc2b')
ax4.plot(times, numeric_4, label='numeric 4', color='#000000')
ax4.grid(color='#e6e8ed', linestyle='-', linewidth=2)
ax4.legend()

# Subplot 5
ax5 = plt.subplot2grid((5, 2), (4,0), colspan=1, rowspan=1, sharex=ax1, sharey=ax1)
ax5.plot(times, analytic_5, label='analytic 5', color='#0306bc')
ax5.plot(times, numeric_5, label='numeric 5', color='#000000')
ax5.grid(color='#e6e8ed', linestyle='-', linewidth=2)
ax5.set_xlabel('Time(hours)')
ax5.legend()

# Subplot 6: the contour plot
ax6 = plt.subplot2grid((5, 2), (0,1), colspan=1, rowspan=5)

# Resolution of the co-latitude
n1 = 1000
# Resolution of the longitude
n2 = 1000

colat = np.linspace(np.pi, 0, n1)
lng = np.linspace(-np.pi, np.pi, n2)
time = 12
X,Y = np.meshgrid(lng, colat)
colat=Y
lng = X
V = Vnz(colat=colat, lng=lng, i=inclination, obq=obliquity, sol_phase=sol_phase, times=time)
I = Inz(colat=colat, lng=lng, i=inclination, obq=obliquity, sol_phase=sol_phase, times=time)
K = V*I
cs = plt.contourf(X,Y,K, cmap='gray', extend='both', vmin=0)
plt.colorbar()

# Converting the obl, xi_sol, and inc values to 2 decimal places for title
obl = "{:.2f}".format(obliquity)
xi_sol = "{:.2f}".format(sol_phase)
inc = "{:.2f}".format(inclination)

plt.title(r'Time: 0{}h, $\Theta$ = {}, $\xi_s$ = {}, $i$ = {}'
          .format(time, obl, xi_sol, inc))
plt.xlabel(r'$\phi$')
plt.ylabel(r'$\theta$')
plt.ylim(np.pi, 0)
plt.xlim(-np.pi, np.pi)

# Setting the ticks to be in LaTeX
plt.xticks([-np.pi, -np.pi/2, 0, np.pi/2, np.pi],
           [r'$- \pi$', r'$- \frac{\pi}{2}$', '$0$',
            r'$\frac{\pi}{2}$', r'$\pi$'])
plt.yticks([np.pi, 3*np.pi/4, np.pi/2, np.pi/4, 0],
           [r'$\pi$', r'$\frac{3 \pi}{4}$', r'$\frac{\pi}{2}$',
            r'$\frac{\pi}{4}$', '$0$'])

# Plotting the points now
plt.scatter([-0.8], [0], c='#f20e0e')
plt.scatter([-0.8], [0.523599], c='#ff7f00')
plt.scatter([-0.8], [0.785398], c='#ffc300')
plt.scatter([-0.8], [1.0472], c='#03bc2b')
plt.scatter([-0.8], [1.5708], c='#0306bc')

# Adding in the contours for I and V separately
CS1 = plt.contour(X,Y,V, 6, colors='r')
CS2 = plt.contour(X,Y,I, 6, colors='b')

fig.tight_layout()
fig.subplots_adjust(hspace=0)

# Saving the figure
title = 'plots/[Side-by-side]CONTOUR PLOT, i={}, obq={}, sph={}, time=0{}h.pdf'.format(inc, obl, xi_sol, time)
# fig.savefig(title)

plt.show()

# Making the contour map only
fig2=plt.figure()
time = time
cs = plt.contourf(X,Y,K, cmap='gray', extend='both', vmin=0)
plt.colorbar()
plt.title(r'Time: 0{}h, $\Theta$ = {}, $\xi_s$ = {}, $i$ = {}'
          .format(time, obl, xi_sol, inc))
plt.xlabel(r'$\phi$')
plt.ylabel(r'$\theta$')
plt.ylim(np.pi, 0)
plt.xlim(-np.pi, np.pi)

# Setting the ticks to be in LaTeX
plt.xticks([-np.pi, -np.pi/2, 0, np.pi/2, np.pi],
           [r'$- \pi$', r'$- \frac{\pi}{2}$', '$0$',
            r'$\frac{\pi}{2}$', r'$\pi$'])
plt.yticks([np.pi, 3*np.pi/4, np.pi/2, np.pi/4, 0],
           [r'$\pi$', r'$\frac{3 \pi}{4}$', r'$\frac{\pi}{2}$',
            r'$\frac{\pi}{4}$', '$0$'])

# Plotting the points now
plt.scatter([-0.8], [0], c='#f20e0e')
plt.scatter([-0.8], [0.523599], c='#ff7f00')
plt.scatter([-0.8], [0.785398], c='#ffc300')
plt.scatter([-0.8], [1.0472], c='#03bc2b')
plt.scatter([-0.8], [1.5708], c='#0306bc')

# Adding in the contours for I and V separately
CS1 = plt.contour(X,Y,V, 6, colors='r')
CS2 = plt.contour(X,Y,I, 6, colors='b')

title = 'plots/[Just contour plot]CONTOUR PLOT, i={}, obq={}, sph={}, time=0{}h.pdf'.format(inc, obl, xi_sol, time)
# fig2.savefig(title)
plt.show()

# Next part of debugging: obtaining phi_obs and phi_stellar. First,
# going to obtain the analytic versions

def get_phi_s(colat, lng, w_rot, w_orb, obq, i, sol_phase, times, nside):
    # (C1)
    cos_theta_knot = (np.cos(i)*np.cos(obq)) + (np.sin(i)*np.sin(obq)*np.cos(sol_phase))
    sin_theta_knot = np.sqrt(1-(cos_theta_knot**2))

    if i == np.pi/2:
        # (C1) Orbital Inclination
        cos_theta_knot = np.sin(obq)*np.cos(sol_phase)

        # (C12)
        cos_phi_s_1 = (np.cos(w_rot*times))*(np.cos(w_orb*times) - cos_theta_knot*np.sin(obq)*np.cos(w_orb*times-sol_phase))
        cos_phi_s_2 = np.sin(w_rot*times)*np.sin(w_orb*times)*np.cos(obq)
        cos_phi_s_3 = sin_theta_knot*(np.sqrt(1-((((np.sin(obq))**2)*((np.cos(w_orb*times-sol_phase))**2)))))
        cos_phi_s = (cos_phi_s_1+cos_phi_s_2)/(cos_phi_s_3)

        # (C13)
        sin_phi_s_1 = -np.sin(w_rot*times)*(np.cos(w_orb*times)-cos_theta_knot*np.sin(obq)*np.cos(w_orb*times-sol_phase))
        sin_phi_s_2 = np.cos(w_rot*times)*np.sin(w_orb*times)*np.cos(obq)
        sin_phi_s_3 = sin_theta_knot*(np.sqrt(1-(((np.sin(obq))**2)*((np.cos(w_orb*times-sol_phase))**2))))
        sin_phi_s = (sin_phi_s_1+sin_phi_s_2)/sin_phi_s_3

    elif i == 0 :

        # (C14)
        cos_phi_s_1 = -np.cos(w_rot*times)*np.cos(obq)*np.cos(w_orb*times-sol_phase)
        cos_phi_s_2 = np.sin(w_rot*times)*np.sin(w_orb*times-sol_phase)
        cos_phi_s_3 = np.sqrt(1-((((np.sin(obq))**2))*((np.cos(w_orb*times-sol_phase))**2)))
        cos_phi_s = (cos_phi_s_1-cos_phi_s_2)/cos_phi_s_3

        # (C15)
        sin_phi_s_1 = np.sin(w_rot*times)*np.cos(obq)*np.cos(w_orb*times-sol_phase)
        sin_phi_s_2 = np.cos(w_rot*times)*np.sin(w_orb*times-sol_phase)
        sin_phi_s_3 = np.sqrt(1-(((np.sin(obq))**2)*((np.cos(w_orb*times-sol_phase))**2)))
        sin_phi_s = (sin_phi_s_1-sin_phi_s_2)/sin_phi_s_3

    elif obq == 0:
        cos_phi_s = np.cos((w_orb-w_rot)*times)
        sin_phi_s = np.sin((w_orb-w_rot)*times)

    elif obq == np.pi/2:
        # (C20)
        cos_theta_knot = np.sin(i)*np.sin(obq)*np.cos(sol_phase)

        # (C22)
        cos_phi_s_1 = np.cos(w_rot*times)*(np.sin(i)*np.cos(w_orb*times) - cos_theta_knot*np.cos(w_orb*times-sol_phase))
        cos_phi_s_2 = np.sin(w_rot*times)*np.cos(i)*np.sin(w_orb*times-sol_phase)
        cos_phi_s_3 = np.sqrt(1-(((np.sin(i))**2)*((np.cos(sol_phase))**2)))
        cos_phi_s_4 = np.sin(w_orb*times-sol_phase)
        cos_phi_s = (cos_phi_s_1-cos_phi_s_2)/(cos_phi_s_3*cos_phi_s_4)

        # (C23)
        sin_phi_s_1 = -np.sin(w_rot*times)*(np.sin(i)*np.cos(w_orb*times) - cos_theta_knot*np.cos(w_orb*times-sol_phase))
        sin_phi_s_2 = np.cos(w_rot*times)*np.cos(i)*np.sin(w_orb*times-sol_phase)
        sin_phi_s_3 = np.sqrt(1-(((np.sin(i))**2)*((np.cos(sol_phase))**2)))
        sin_phi_s_4 = np.sin(w_orb*times-sol_phase)
        sin_phi_s = (sin_phi_s_1-sin_phi_s_2)/(sin_phi_s_3*sin_phi_s_4)

    else:
        # (C4)
        cos_phi_s_1 = np.cos(w_rot*times)*(np.sin(i)*np.cos(w_orb*times) - cos_theta_knot*np.sin(obq)*np.cos(w_orb*times-sol_phase))
        cos_phi_s_2 = np.sin(w_rot*times)*(np.sin(i)*np.sin(w_orb*times)*np.cos(obq) - np.cos(i)*np.sin(obq)*np.sin(w_orb*times-sol_phase))
        cos_phi_s_3 = np.sqrt(1-(np.cos(i)*np.cos(obq) + np.sin(i)*np.sin(obq)*np.cos(sol_phase))**2)
        cos_phi_s_4 = np.sqrt(1-(((np.sin(obq))**2)*((np.cos(w_orb*times-sol_phase))**2)))
        cos_phi_s = (cos_phi_s_1+cos_phi_s_2)/(cos_phi_s_3*cos_phi_s_4)

        # (C5)
        sin_phi_s_1 = -np.sin(w_rot*times)*(np.sin(i)*np.cos(w_orb*times) - cos_theta_knot*np.sin(obq)*np.cos(w_orb*times-sol_phase))
        sin_phi_s_2 = np.cos(w_rot*times)*(np.sin(i)*np.sin(w_orb*times)*np.cos(obq) - np.cos(i)*np.sin(obq)*np.sin(w_orb*times-sol_phase))
        sin_phi_s_3 = np.sqrt(1-(np.cos(i)*np.cos(obq)+np.sin(i)*np.sin(obq)*np.cos(sol_phase))**2)
        sin_phi_s_4 = np.sqrt(1-(((np.sin(obq))**2)*((np.cos(w_orb*times-sol_phase))**2)))
        sin_phi_s = (sin_phi_s_1+sin_phi_s_2)/(sin_phi_s_3*sin_phi_s_4)
    return sin_phi_s, cos_phi_s

def get_phi_o(colat, lng, w_rot, w_orb, obq, i, sol_phase, times, nside):
    cos_phi_knot = np.cos(w_rot*times)
    sin_phi_knot = -np.sin(w_rot*times)
    return sin_phi_knot, cos_phi_knot


sin_phi_s,cos_phi_s = get_phi_s(colat=colat_1, lng=lng, w_rot=w_rot,
                                w_orb=w_orb, obq=obliquity, i=inclination, sol_phase=sol_phase,
                                times=time, nside=nside)
sin_phi_knot, cos_phi_knot = get_phi_o(colat=colat_1, lng=lng, w_rot = w_rot,
                                       w_orb=w_orb, obq=obliquity, i=inclination, sol_phase=sol_phase,
                                       times=time, nside=nside)

phi_s, phi_knot = 0, 0

# Using try-catch statements to obtain the values for phi_s and phi_knot
try:
    phi_s = np.arcsin(sin_phi_s)
except:
    phi_s = np.arccos(cos_phi_s)

try:
    phi_knot = np.arcsin(sin_phi_knot)
except:
    phi_knot = np.arccos(cos_phi_knot)


# Now, obtaining the trig values for the numeric solution
xi_0 = 0
time = np.array([time])
trigvals = viewgeom(times=time, wrot=w_rot, worb=w_orb, obq=obliquity, inc=inclination,
                    xisol=sol_phase, xi0=xi_0)

n_sin_phi_s = trigvals.item(6)
n_cos_phi_s = trigvals.item(7)
n_sin_phi_o = trigvals.item(2)
n_cos_phi_o = trigvals.item(3)


# Obtaining the values for phi_s and phi_o
n_phi_s, n_phi_o = 0, 0

try:
    n_phi_s = np.arcsin(n_sin_phi_s)
except:
    n_phi_s = np.arccos(n_cos_phi_s)

try:
    n_phi_o = np.arcsin(n_sin_phi_o)
except:
    n_phi_o = np.arccos(n_cos_phi_o)

# Now going to obtain theta_s and theta_o since there is perfet agreement
# with the phi_s and phi_o

def get_theta_s(colat, lng, w_rot, w_orb, obq, i, sol_phase, times, nside):
    cos_theta_s = np.sin(obq)*np.cos((w_orb*times)-sol_phase)
    sin_theta_s = np.sqrt(1-(((np.sin(obq))**2)*((np.cos((w_orb*times)-sol_phase))**2)))
    return sin_theta_s, cos_theta_s

def get_theta_o(colat, lng, w_rot, w_orb, obq, i, sol_phase, times, nside):
    cos_theta_knot = (np.cos(i)*np.cos(obq)) + (np.sin(i)*np.sin(obq)*np.cos(sol_phase))
    sin_theta_knot = np.sqrt(1-(cos_theta_knot**2))
    return sin_theta_knot, cos_theta_knot


sin_theta_s, cos_theta_s = get_theta_s(colat=colat_1, lng=lng, w_rot=w_rot,
                                       w_orb=w_orb, obq=obliquity,
                                       i=inclination, sol_phase=sol_phase,
                                       times=time, nside=nside)
sin_theta_o, cos_theta_o = get_theta_o(colat=colat_1, lng=lng, w_rot=w_rot,
                                       w_orb=w_orb, obq=obliquity, i=inclination,
                                       sol_phase=sol_phase, times=time, nside=nside)

# Obtaining the values for theta_s and theta_o
theta_s, theta_o = 0, 0

try:
    theta_s = np.arcsin(sin_theta_s)
except:
    theta_s = np.arccos(cos_theta_s)

try:
    theta_o = np.arcsin(sin_theta_o)
except:
    theta_o = np.arccos(cos_theta_o)

# Obtaining the numeric solution's values for theta_s and theta_o
n_sin_theta_o = trigvals.item(0)
n_cos_theta_o = trigvals.item(1)
n_sin_theta_s = trigvals.item(4)
n_cos_theta_s = trigvals.item(5)

# Obtaining the values for theta_s and theta_o

n_theta_o, n_theta_s = 0,0
try:
    n_theta_o = np.arcsin(sin_theta_o)
except:
    n_theta_o = np.arccos(cos_theta_o)
try:
    n_theta_s = np.arcsin(sin_theta_s)
except:
    n_theta_s = np.arccos(cos_theta_s)



# # Printing the results
# print("Trig value comparison results:")
# print("cos_theta_s [numeric]: {} | [analytic]: {}".format(n_cos_theta_s, cos_theta_s))
# print("sin_theta_s [numeric]: {} | [analytic]: {}".format(n_sin_theta_s, sin_theta_s))
# print("theta_s [numeric]: {}     | [analytic]: {}".format(n_theta_s, theta_s))
# print("sin_theta_o [numeric]: {} | [analytic]: {}".format(n_sin_theta_o, sin_theta_o))
# print("cos_theta_o [numeric]: {} | [analytic]: {}".format(n_cos_theta_o, cos_theta_o))
# print("theta_o [numeric]: {}     | [analytic]: {}".format(n_theta_o, theta_o))
# print("cos_phi_s: [numeric]: {}  | [analytic]: {}".format(n_cos_phi_s, cos_phi_s))
# print("sin_phi_s: [numeric]: {}  | [analytic]: {}".format(n_sin_phi_s, sin_phi_s))
# print("phi_s: [numeric]: {}      | [analytic]: {}".format(n_phi_s, phi_s))
# print("cos_phi_o: [numeric]: {}  | [analytic]: {}".format(n_cos_phi_o, cos_phi_knot))
# print("sin_phi_o: [numeric]: {}  | [analytic]: {}".format(n_sin_phi_o, sin_phi_knot))

is_close_cos_theta_s = math.isclose(n_cos_theta_s, cos_theta_s, rel_tol=1e-5)
is_close_sin_theta_s = math.isclose(n_sin_theta_s, sin_theta_s, rel_tol=1e-5)
is_close_theta_s = math.isclose(n_theta_s, theta_s, rel_tol=1e-5)

is_close_sin_theta_o = math.isclose(n_sin_theta_o, sin_theta_o, rel_tol=1e-5)
is_close_cos_theta_o = math.isclose(n_cos_theta_o, cos_theta_o, rel_tol=1e-5)
is_close_theta_o = math.isclose(n_theta_o, theta_o, rel_tol=1e-5)

is_close_cos_phi_s = math.isclose(n_cos_phi_s, cos_phi_s, rel_tol=1e-5)
is_close_sin_phi_s = math.isclose(n_sin_phi_s, sin_phi_s, rel_tol=1e-5)
is_close_phi_s = math.isclose(n_phi_s, phi_s, rel_tol=1e-5)

is_close_cos_phi_o = math.isclose(n_cos_phi_o, cos_phi_knot, rel_tol=1e-5)
is_close_sin_phi_o = math.isclose(n_sin_phi_o, sin_phi_knot, rel_tol=1e-5)
is_close_phi_o = math.isclose(n_phi_o, phi_knot, rel_tol=1e-5)

boolean= is_close_theta_s and is_close_theta_o and is_close_phi_s and is_close_phi_o
print("Do all the angles agree? {}".format(boolean))

# Now, going to produce a contour map for the *NUMERIC* solution

phi_rot = 0
times = 12
V, I = get_visibility_illum(p_rotation=p_rotation, p_orbit=p_orbit, phi_orb=phi_orb,
                            inclination=inclination, solstice_phase=sol_phase,
                            obliquity=obliquity, phi_rot =phi_rot, times=times)
K = V*I

# Obtaining the sub-stellar and sub-observer points
nos, nss = get_points(p_rotation=p_rotation, p_orbit=p_orbit, phi_orb=phi_orb,
                      inclination=inclination, solstice_phase=sol_phase, obliquity=obliquity,
                      phi_rot=phi_rot, times=times)

# Now converting the cartesian coordinates to spherical coordinates

def cartesian_to_spherical(cartesian):
    x = cartesian.item(0)
    y = cartesian.item(1)
    z = cartesian.item(2)

    r = np.sqrt(x*x + y*y + z*z)
    theta = np.arccos(z/r)
    phi = np.arctan(y/x)
    return theta, phi

nos_spherical_theta, nos_spherical_phi = cartesian_to_spherical(nos)
nss_spherical_theta, nss_spherical_phi = cartesian_to_spherical(nss)


# Printing everything together:
print("Sub-stellar point: ")
print("Theta_s: [Analytic]: {} [Numeric]: {}".format(theta_s, nss_spherical_theta))
print("Phi_s: [Analytic]: {} [Numeric]: {}".format(phi_s, nss_spherical_phi))
print("Sub-observer point: ")
print("Theta_o: [Analytic]: {} [Numeric]: {}".format(theta_o, nos_spherical_theta))
print("Phi_o: [Analytic]: {} [Numeric]: {}".format(phi_knot, nos_spherical_phi))