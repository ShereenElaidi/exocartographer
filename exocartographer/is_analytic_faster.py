import time as t
from analytic_delta import analytic_delta
import numpy as np
import matplotlib.pyplot as plt
from sklearn.metrics import mean_squared_error
from scipy.spatial import distance
from math import sqrt
from numerically_integrate import numerically_integrate
import healpy as hp

# SIMULATING AN ENTIRE PLANET'S LIGHTCURVE
def make_lightcurve_plot(flux):
    plt.plot(flux)
    plt.ylabel('flux')
    plt.xlabel('time(rotations)')
    plt.show()

# This will simulate an entire planet
theta = [0,np.pi/6, np.pi/5, np.pi/4, np.pi/3, np.pi/2,
         2*np.pi/3, 3*np.pi/4, 5*np.pi/6, np.pi]
phi = [0,np.pi/6, np.pi/5, np.pi/4, np.pi/3, np.pi/2,
         2*np.pi/3, 3*np.pi/4, 5*np.pi/6, np.pi]


time_data = np.array([0, 1, 2, 3, 4, 5, 6, 7, 8, 9, 10])
w_rot = 1
w_orb = 1
obliq = 1
i = 1
sol_phase = 1
orb_pos = 1
sub_long = 1
colat = np.pi/3                    # location of bright point
colng = np.pi/3                    # location of bright point
size = 2                            # pixel size for numerical integration
xi0 = 1

# -----------------CALCULATING THE TIME DIFFERENCE---------------------------------

# Put code for analytic in here:
t0 = t.clock()
for k in range(0, len(theta)):
    for j in range(0, len(phi)):
        flux = analytic_delta(w_rot, w_orb, obliq, i, sol_phase, orb_pos,
                              sub_long, time_data, int(theta[k]), int(phi[j]))
t1 = t.clock()
analytic_time = t1-t0

# Put code for numerical integration here:

# def numerically_integrate(lat, lng, times, wrot, worb, obq, inc, xisol, xi0, colat, colng, size):
t2 = t.clock()
for k in range(0, len(theta)):
    for j in range(0, len(phi)):
        flux = numerically_integrate(theta[k], phi[j], time_data, w_rot, w_orb, obliq, i, sol_phase, xi0,  colat, colng, size)

t3 = t.clock()
numerical_time = t3 - t2

print("Computing values....")

# -------CALCULATING THE ACCURACY OF THE ANALYTIC MAP-----------------

# True values - the values of the flux for the numerical integration:
true_flux = []

for k in range(0, len(theta)):
    for j in range(0, len(phi)):
        true_flux.append(numerically_integrate(theta[k], phi[j], time_data, w_rot, w_orb, obliq, i, sol_phase, orb_pos,
                                               colat, colng, size))

# Estimated values - the estimated flux values w/ analytic
analytic_flux = []
for k in range(0, len(theta)):
    for j in range(0, len(phi)):
        analytic_flux.append(analytic_delta(w_rot, w_orb, obliq, i, sol_phase, orb_pos,
                              sub_long, time_data, theta[k], phi[j]))


def calculate_rmse(true_flux, analytic_flux):
    return sqrt(mean_squared_error(true_flux, analytic_flux))


def euclidean_distance(true_flux, analytic_flux, times):
    euclidean_distance = 0
    for row in range(0, len(true_flux)):
        euclidean_distance+= distance.euclidean(true_flux[row], analytic_flux[row])
    return euclidean_distance, euclidean_distance/(len(times)*len(true_flux))


print("-------------------SUMMARY------------------------")
print('TIME STATISTICS:')
print("Analytic time: %f (seconds)" %(analytic_time))
print("Numeric integration time: %f (seconds)" % (numerical_time))
print("ACCURACY STATISTICS")
print("RMSE: %f" %(calculate_rmse(true_flux, analytic_flux)))
euclideandistance, average_euclidean_distance = euclidean_distance(true_flux, analytic_flux, time_data)
print("Euclidean distance: %f \nAverage Euclidean Distance: %f" % (euclideandistance, average_euclidean_distance))
print("Making plots...")
make_lightcurve_plot(true_flux)
make_lightcurve_plot(analytic_flux)
