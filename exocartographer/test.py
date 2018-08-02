import numpy as np
from analytic_lc import analytic_lightcurve
import matplotlib.pyplot as plt

# def analytic_lightcurve(lat, lng, w_rot, w_orb, obq, i, sol_plase, times):

lat = np.pi/4
lng = np.pi/3
w_rot= (2*np.pi)/1
w_orb = 2*np.pi/365
obq = np.pi/4
i = 1
sol_phase = 1
times = np.arange(0.0, 10.0, 0.01)

light_curve = analytic_lightcurve(lat=lat, lng=lng, w_rot=w_rot, w_orb=w_orb,obq= obq, i=i, sol_phase=sol_phase, times=times)
# print(light_curve)

# plt.subplot(2, 1, 2)
# # plt.figure(figsize=(16, 3))
# plt.title("Analytic lightcurve")
# plt.plot(time_data,flux)
# plt.ylabel('flux')
# plt.xlabel(r'time$/P_\mathrm{rot}$')
# plt.show()

plt.title("Analytic lightcurve")
plt.plot(times,light_curve)
plt.ylabel('reflectance')
plt.xlabel('time(rotations)')
plt.show()