# Import statemennts
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

# Broad overview od what the code is going to do :
# (1) obtain the angle values from the *analytic* VIK from
# exocartographer
# (2) re-obtain the angle values from the analytic solution's
# observational angles
# (3) re-obtain the angle values from the numeric solution's angles
# in exocartographer
# (4) plot all these together. I think that the way that I was
# retrieving some of the angles caused some of them to be incorrect

# Helper functions : a function to draw the HEALPIx map
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

# Setting the orbital properties for the numeric solution
nside = 4

# From the numeric solution : exocartographer uses order = ring
order = 'ring'
map = HEALPix(nside=nside, order=order)

# The unit here will be in h
p_rotation = 23.934
p_orbit = 365.256363 * 24.0

# Creating an evenly-spaced time array
times = np.linspace(start=0.0, stop=365*24.0, num=14000)
measurement_std = 0.001

# Setting the input parameters :
w_rot = 2 * np.pi/p_rotation
w_orb = 2 * np.pi/p_orbit
inclination = np.pi/4
obliquity = np.pi/6
sol_phase = np.pi/6

# Using the relation between exocartographer's inputs and the numeric
# code's inputs here
phi_orb = 0                         # xi_0
obl_orientation = sol_phase         # xi_s

# Creating the point of the ``bright spot''
lat_1 = Angle((np.pi/2) * u.rad)
lng = Angle((-np.pi/2)*u.rad)

# Converting the latitude value into co-latitude values
a = Angle((np.pi/2) * u.rad)
colat_1 = a-lat_1

# Obtaining the HEALPix index for the bright spot
bs = map.lonlat_to_healpix([lng]*u.rad, [lat_1]*u.rad)

# Drawing the simulated map for the bright spot
sim_map = draw_albedo_map(nside=nside, view_map=False, bright_spot=bs)


# Defining all the angles that we need :

# The vector-rotating EXOCARTOGRAPHER's code
e_phi_o = np.zeros_like(times)
e_phi_s = np.zeros_like(times)
e_th_o = np.zeros_like(times)
e_th_s = np.zeros_like(times)

# My analytic solution's code
a_phi_o = np.zeros_like(times)
a_phi_s = np.zeros_like(times)
a_th_o = np.zeros_like(times)
a_th_s = np.zeros_like(times)

# The analytic code from exocartographer
ea_phi_o = np.zeros_like(times)
ea_phi_s = np.zeros_like(times)
ea_th_o = np.zeros_like(times)
ea_th_s = np.zeros_like(times)

# ANALYIC ANGLES

# You say phi_s should change linearly with time, but the way that
# this angle is retrieved means that we're effectively mod(2pi)-ing the
# angle's value...not sure how to counter-act this since there's no explicit
# formula for phi_s, unlike phi_o

def get_phi_s( w_rot, w_orb, obq, i, sol_phase, times):
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
    return np.arcsin(sin_phi_s)

# Same issue as earlier with this one, so I just re-wrote the function
# so that it does not get modded.
def get_phi_o(w_rot, times):
    s_phi_o = -np.sin(w_rot*times)
    phi_o = np.arcsin(s_phi_o)
    return phi_o

def get_theta_s(obq, w_orb, times, sol_phase):
    sin_theta_s = np.sqrt(1 - (((np.sin(obq)) ** 2) * ((np.cos((w_orb * times) - sol_phase)) ** 2)))
    return np.arcsin(sin_theta_s)

def get_theta_o(i, obq, sol_phase, times):
    cos_theta_knot = (np.cos(i)*np.cos(obq)) + (np.sin(i)*np.sin(obq)*np.cos(sol_phase))
    return np.full_like(times, np.arccos(cos_theta_knot))

a_phi_s = get_phi_s(w_rot=w_rot, w_orb=w_orb, obq=obliquity,
                              i=inclination, sol_phase=sol_phase,
                              times=times)
a_phi_o = get_phi_o(w_rot=w_rot, times=times)

a_th_s = get_theta_s(obq=obliquity, w_orb=w_orb, times=times,
                     sol_phase=sol_phase)
a_th_o = get_theta_o(i=inclination, obq=obliquity, sol_phase=sol_phase, times=times)


# NUMERIC SOLUTION'S VECTOR ROTATION CODE :
xi_0 = 0
trigvals = viewgeom(times=times, wrot=w_rot, worb=w_orb, obq=obliquity,
                   inc=inclination, xisol=sol_phase, xi0=xi_0)

ea_phi_s = np.arcsin(trigvals[:,6]) #6th
ea_phi_o = np.arcsin(trigvals[:,2]) #2nd
ea_th_s = np.arcsin(trigvals[:,4])  #4th
ea_th_o = np.arcsin(trigvals[:,0])  #0th


# EXOCARTOGRAPHER'S NUMERIC VERSION :
nos, nss = get_points(p_rotation = p_rotation, p_orbit=p_orbit,
                      phi_orb=phi_orb, inclination=inclination,
                      solstice_phase=sol_phase, obliquity=obliquity,
                      phi_rot=phi_orb, times=times)

# Extracting the th_o angle
c_th_o = nos[:,2]
e_th_o = np.arccos(c_th_o)

# Extracting the th_s angle
c_th_s = nss[:,2]
e_th_s = np.arccos(c_th_s)


# Extracting the phi_o angle
s_th_o = np.sin(e_th_o)
e_phi_o = np.arcsin(nos[:,1]/s_th_o)

# Extracting the phi_s angle
s_th_s = np.sin(e_th_s)
e_phi_s = np.arcsin(nss[:,1]/s_th_s)

# Plotting everything together

f, (ax1, ax2, ax3, ax4) = plt.subplots(4, 1, sharex=True, sharey=True)
ax1.plot(times, ea_phi_s, label='analytic - exocartographer')
ax1.plot(times, a_phi_s, label='analytic - my code')
ax1.plot(times, e_phi_s, label='numeric')
ax1.grid(color='#e6e8ed', linestyle='-', linewidth=2)
ax1.legend()
ax1.set_ylabel(r'$\phi_s$')

ax2.plot(times, ea_phi_o, label='analytic - exocartographer')
ax2.plot(times, a_phi_o, label='analytic - my code')
ax2.plot(times, e_phi_o, label='numeric')
ax2.grid(color='#e6e8ed', linestyle='-', linewidth=2)
ax2.legend()
ax2.set_ylabel(r'$\phi_o$')

ax3.plot(times, ea_th_s, label='analytic - exocartographer')
ax3.plot(times, a_th_s, label='analytic - my code')
ax3.plot(times, e_th_s, label='numeric')
ax3.grid(color='#e6e8ed', linestyle='-', linewidth=2)
ax3.legend()
ax3.set_ylabel(r'$\theta_s$')

ax4.plot(times, ea_th_o, label='analytic - exocartographer')
ax4.plot(times, a_th_o, label='analytic - my code')
ax4.plot(times, e_th_o, label='numeric')
ax4.grid(color='#e6e8ed', linestyle='-', linewidth=2)
ax4.legend()
ax4.set_ylabel(r'$\theta_o$')
ax4.set_xlabel('time(h)')

f.tight_layout()
f.subplots_adjust(hspace=0)
plt.show()