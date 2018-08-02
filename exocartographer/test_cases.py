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
import pdb

# A functon to plot curves
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


# Drawing the map
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

# Setting the orbital properties
nside = 4
order = 'ring'
map = HEALPix(nside=nside, order=order)

# The unit here will be hours
p_rotation = 23.934
p_orbit = 365.256363 * 24.0

# The unit here will be rotations
# p_rotation = 1
# p_orbit =  365.256363

times = np.linspace(start=0.0, stop=24.0, num=1400)
measurement_std = 0.001
w_rot = 2*np.pi/p_rotation
w_orb = 2*np.pi/p_orbit
inclination = 0
obliquity = 0
sol_phase = 0  # in [0, 2pi[

# Relation (1)
# phi_orb = abs((7*(np.pi/6)) + sol_phase)
# phi_rot = abs((np.pi/3) + sol_phase)

# Relation (2)
# phi_orb = abs((2*(np.pi/3)) + sol_phase)
# phi_rot = abs((2*(np.pi/3)) + sol_phase)

# Relation (3)
phi_orb = abs((1*(np.pi/4)) - sol_phase)
phi_rot = abs((1*(np.pi/4)) + sol_phase)
# p_rotation = 2*np.pi/w_rot
# p_orbit = 2*np.pi/w_orb

# Drawing the simulated map. Same longitude, but varying co-latitudes
print("Drawing the simulated map. Same longitude, but varying co-latitudes...")
lat_1 = Angle((np.pi/2)*u.rad)
lat_2 = Angle((np.pi/3)*u.rad)
lat_3 = Angle((np.pi/4)*u.rad)
lat_4 = Angle((np.pi/6)*u.rad)
lat_5 = Angle((0)*u.rad)

# Converting the latitude values into co-latitude values
print("Converting the latitude values into co-latitude values...")
a = Angle((np.pi/2)*u.rad)
colat_1 = a-lat_1
colat_2 = a-lat_2
colat_3 = a-lat_3
colat_4 = a-lat_4
colat_5 = a-lat_5

lng = Angle(-(np.pi/2)*u.rad)

# Obtaining the bright spot for each point
print("Obtaining the bright spot for each point...")
print("Bright point 1...")
bs_1 = map.lonlat_to_healpix([lng]*u.rad, [lat_1]*u.rad)
print("Bright point 2...")
bs_2 = map.lonlat_to_healpix([lng]*u.rad, [lat_2]*u.rad)
print("Bright point 3...")
bs_3 = map.lonlat_to_healpix([lng]*u.rad, [lat_3]*u.rad)
print("Bright point 4...")
bs_4 = map.lonlat_to_healpix([lng]*u.rad, [lat_4]*u.rad)
print("Bright point 5...")
bs_5 = map.lonlat_to_healpix([lng]*u.rad, [lat_5]*u.rad)

# Drawing the simulated map for each bright spot
print("Drawing the simulated map for each bright spot")
print("Drawing simulated map 1...")
sim_map_1 = draw_albedo_map(nside=nside, view_map=False, bright_spot=bs_1)
print("Drawing simulated map 2...")
sim_map_2 = draw_albedo_map(nside=nside, view_map=False, bright_spot=bs_2)
print("Drawing simulated map 3...")
sim_map_3 = draw_albedo_map(nside=nside, view_map=False, bright_spot=bs_3)
print("Drawing simulated map 4...")
sim_map_4 = draw_albedo_map(nside=nside, view_map=False, bright_spot=bs_4)
print("Drawing simulated map 5...")
sim_map_5 = draw_albedo_map(nside=nside, view_map=False, bright_spot=bs_5)

# Analytic solution
analytic_1 = analytic_lightcurve_2(colat=colat_1, lng=lng, w_rot=w_rot, w_orb = w_orb,
                                   obq = obliquity, i=inclination, sol_phase=sol_phase,
                                   times=times, nside=nside)
analytic_2 = analytic_lightcurve_2(colat=colat_2, lng=lng, w_rot=w_rot, w_orb = w_orb,
                                   obq=obliquity, i=inclination, sol_phase=sol_phase,
                                   times=times, nside=nside)
analytic_3 = analytic_lightcurve_2(colat=colat_3, lng=lng, w_rot =w_rot, w_orb=w_orb,
                                   obq=obliquity, i=inclination, sol_phase=sol_phase,
                                   times=times, nside=nside)
analytic_4 = analytic_lightcurve_2(colat=colat_4, lng=lng, w_rot=w_rot, w_orb=w_orb,
                                   obq=obliquity, i=inclination, sol_phase=sol_phase,
                                   times=times, nside=nside)
analytic_5 = analytic_lightcurve_2(colat=colat_5, lng=lng, w_rot=w_rot, w_orb=w_orb,
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
    'logit_obl_orientation':logit(phi_rot, low=0, high=2*np.pi)
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

# Plotting all the curves together
print("making plots...")
plt.subplots(1,2)

# Subplot 1
plt.subplot(1,2,1)
f, (ax1, ax2, ax3, ax4, ax5) = plt.subplots(5, sharex=True, sharey=True)
f.suptitle('Analytic and numeric lightcurves at different locations')
ax1.plot(times, analytic_1, label='analytic 1')
ax1.plot(times, numeric_1, label='numeric 1')
ax1.legend()
ax1.set_ylabel('reflectance')
ax1.set_xlabel('time(hours)')

ax2.plot(times, analytic_2, label='analytic 2')
ax2.plot(times, numeric_2, label='numeric 2')
ax2.legend()
ax2.set_ylabel('reflectance')
ax2.set_xlabel('time(hours)')

ax3.plot(times, analytic_3, label='analytic 3')
ax3.plot(times, numeric_3, label='numeric 3')
ax3.legend()
ax3.set_ylabel('reflectance')
ax3.set_xlabel('time(hours)')

ax4.plot(times, analytic_4, label='analytic 4')
ax4.plot(times, numeric_4, label='numeric 4')
ax4.legend()
ax4.set_ylabel('reflectance')
ax4.set_xlabel('time(hours)')

ax5.plot(times, analytic_5, label='analytic 5')
ax5.plot(times, numeric_5, label='numeric 5')
ax5.legend()
ax5.set_ylabel('reflectance')
ax5.set_xlabel('time(hours)')

f.subplots_adjust(hspace=0)
plt.show()

# Making the contour plot for a 24-hour period
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

# Subplot 2

plt.subplot(1,2,2)
# Resolution of co-latitude
n1 = 1000

# Resolution of the longitude
n2 = 1000

colat = np.linspace(np.pi, 0, n1)
lng = np.linspace(-np.pi, np.pi, n2)

# Setting the parameters
i = inclination
obq = obliquity
sol_phase = sol_phase
time = 15.0

# Making the contour plot (bright spot shows the kernel).
# Also sub-plot 2

times=15
X,Y = np.meshgrid(lng, colat)
colat = X
lng = Y
V = Vnz(colat=colat, lng=lng, i=i, obq=obq, sol_phase=sol_phase, times=times)
I = Inz(colat=colat, lng=lng, i=i, obq=obq, sol_phase=sol_phase, times=times)
K = V*I
cs6=plt.contourf(X,Y,K, cmap = 'gray', extend='both', vmin=0)
plt.colorbar()
plt.title('Time: 15h')
plt.xlabel('longitude')
plt.ylabel('colatitude')


plt.show()


