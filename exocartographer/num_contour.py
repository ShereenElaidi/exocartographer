import numpy as np
import healpy as hp
from matplotlib import pyplot as plt
from exocartographer.gp_map import draw_map
from exocartographer import IlluminationMapPosterior
from exocartographer.util import logit, inv_logit
import exocartographer
import pdb
from exocartographer import viewgeom, kernel
import time
# Map resolution
nside = 4

# Gaussian Process Properties
whitenoise_relative_amp = 0.02
length_scale = 30.*np.pi/180
albedo_mean = 0.5
albedo_std = 0.2

# Drawing a valid albedo map
while True:
    simulated_map = draw_map(nside, albedo_mean, albedo_std,
                             whitenoise_relative_amp, length_scale)
    if min(simulated_map) > 0 and max(simulated_map) < 1:
        break

hp.mollview(simulated_map, title='albedo', cmap='gist_gray')

# setting the orbital properties
p_rotation = 23.934
p_orbit = 365.256363 * 24.0
phi_orb = np.pi
phi_rot = np.pi
inclination = 0
obliquity = 0

times = 7

# Measurement uncertainty
measurement_std = 0.001

# Using a posterior instance to generate a light-curve
truth = IlluminationMapPosterior(times, np.zeros_like(times),
                                 measurement_std, nside=nside)
true_params = {
    'log_orbital_period': np.log(p_orbit),
    'log_rotation_period': np.log(p_rotation),
    'logit_cos_inc': logit(np.cos(inclination)),
    'logit_cos_obl': logit(np.cos(obliquity)),
    'logit_phi_orb': logit(phi_orb, low=0, high=2*np.pi),
    'logit_obl_orientation': logit(phi_rot, low=0, high=2*np.pi)
}

# truth.fix_params(true_params)
# p = np.concatenate([np.zeros(truth.nparams), simulated_map])
# kernel = truth.visibility_illumination_matrix(p)
# print(simulated_map)


# Now, trying to make the contour plot
# plt.contourf(kernel)
# plt.show()




# NO: something else. Try the kernel. Code lifted from:
# https://github.com/bfarr/exocartographer/blob/master/exocartographer/analytic_kernel.py

lat = np.pi/2
lng = np.pi/2
times  = np.array([7])
times = np.linspace(0.0, 10.0, 1400)


theta = lat
phi = lng
omega_rot = 2.0*np.pi/p_rotation
obl = obliquity
obl_orientation = phi_rot # ???
omega_orb = 2.0*np.pi/p_orbit
phi_orb = phi_orb
inc = inclination


def numerical_kernel(colat, lng, w_rot, obl, obl_orientation, w_orb, inc, times, phi_orb):
    theta = colat
    phi = lng
    omega_rot = w_rot
    obl = obl
    obl_orientation = obl_orientation
    omega_orb = w_orb
    phi_orb = phi_orb
    inc = inc

    times = times

    # Obtaining the sub_obs and sub_st trig values
    all_trigs = viewgeom(times, omega_rot, omega_orb, obl, inc, obl_orientation, phi_orb)

    sinth = np.sin(theta)
    costh = np.cos(theta)
    sinph = np.sin(phi)
    cosph = np.cos(phi)

    num_kernel = kernel(sinth, costh, sinph, cosph, all_trigs).T
    return num_kernel

# Computing running-time
t0 = time.clock()

n1 = 100
n2 = 100
colat = np.linspace(np.pi, 0, n1)
lng = np.linspace(-np.pi, np.pi, n2)
X,Y = np.meshgrid(lng, colat)
colat = X
lng = Y
w_rot = 2.0*np.pi/p_rotation
w_orb = 2.0*np.pi/p_orbit
sol_phase = np.pi/6
phi_orb = sol_phase
phi_rot = sol_phase
times = np.array([7])

Z = numerical_kernel(colat=X, lng = Y, w_rot=w_rot, obl = obl, obl_orientation=obl_orientation,
                     w_orb=w_orb, inc = inc, times=times, phi_orb = phi_orb)
print(Z.size)
t1 = time.clock()
time = t1-t0
print("run-time: {}".format(time))

times = 7

cs1 = plt.contourf(X,Y,Z, cmap='gray', extend='both', vmin=0)
plt.colorbar()
plt.title('Time: 0{}'.format(times))
plt.xlabel('longitude')
plt.ylabel('colatitude')