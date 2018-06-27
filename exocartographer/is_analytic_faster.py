import time as t
from analytic_delta import analytic_delta
import numpy as np
import matplotlib.pyplot as plt
from sklearn.metrics import mean_squared_error
from math import sqrt
from numerically_integrate import numerically_integrate

# SIMULATING AN ENTIRE PLANET'S LIGHTCURVE
def make_lightcurve_plot(flux):
    plt.plot(flux)
    plt.ylabel('flux')
    plt.xlabel('time(rotations)')
    plt.show()

# Array of Lat (theta) and Lng (phi)
theta = [0,np.pi/6, np.pi/5, np.pi/4, np.pi/3, np.pi/2,
         2*np.pi/3, 3*np.pi/4, 5*np.pi/6, np.pi, 7*np.pi/6,
         5*np.pi/4,4*np.pi/3, 5*np.pi/3, 7*np.pi/4, 11*np.pi/6]
phi = [0,np.pi/6, np.pi/5, np.pi/4, np.pi/3, np.pi/2,
         2*np.pi/3, 3*np.pi/4, 5*np.pi/6, np.pi, 7*np.pi/6,
         5*np.pi/4,4*np.pi/3, 5*np.pi/3, 7*np.pi/4, 11*np.pi/6]


time_data = np.array([0,1,2,3,4,5,6,7,8,9,10])
w_rot = 1
w_orb = 1
obliq = 1
i = 1
sol_phase = 1
orb_pos = 1
sub_long = 1


# CALCULATING THE TIME DIFFERENCE

# Put code for analytic in here
t0 = t.clock()
for k in range(0, len(theta)):
    for j in range(0, len(phi)):
        flux = analytic_delta(w_rot, w_orb, obliq, i, sol_phase, orb_pos,
                              sub_long, time_data, int(theta[k]), int(phi[j]))
t1 = t.clock()
analytic_time = t1-t0

# Put code for numerical integration here
t2 = t.clock()
for k in range(0, len(theta)):
    for j in range(0, len(phi)):
        flux = analytic_delta(w_rot, w_orb, obliq, i, sol_phase, orb_pos,
                              sub_long, time_data, int(theta[k]), int(phi[j]))
t3 = t.clock()
numerical_time = t3 - t2

print("Analytic time: %f" %(analytic_time))
print("Numeric integration time: %f" % (numerical_time))

# CALCULATING THE ACCURACY OF THE ANALYTIC MAP

# True values - the values of the flux for the numerical integration:
true_flux = []

# Filler code for now until talk to Nick



for k in range(0, len(theta)):
    for j in range(0, len(phi)):
        true_flux.append(numerically_integrate(theta[i], phi[j], time_data, w_rot, w_orb, obliq, i, sol_phase, orb_pos))

# Estimated values - the estimated flux values w/ analytic
analytic_flux = []
for k in range(0, len(theta)):
    for j in range(0, len(phi)):
        analytic_flux.append(analytic_delta(w_rot, w_orb, obliq, i, sol_phase, orb_pos,
                              sub_long, time_data, int(theta[k]), int(phi[j])))


def calculate_rmse(true_flux, analytic_flux):
    return sqrt(mean_squared_error(true_flux, analytic_flux))


print("-------------------SUMMARY------------------------")
print('Statistics for time: ')
print("Analytic time: %f" %(analytic_time))
print("Numeric integration time: %f" % (numerical_time))
print("Accuracy statistics:")
print("RMSE: %f" %(calculate_rmse(true_flux, analytic_flux)))