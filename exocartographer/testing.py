import matplotlib.pyplot as plt
from analytic_delta import analytic_delta
import numpy as np

time = []
y = 0
for y in range (0, 10):
    y = y/10.0
    time.append(y)
    print(y)
print(time)

w_rot = 1
w_orb = 1
obliq = 1
i =1
sol_phase = 1
orb_pos = 1
sub_long = 1
theta_1 = np.pi/4
phi = np.pi/4

theta_2 = np.pi/5
theta_3 = np.pi/5
theta_4 = np.pi/6

array_1 = analytic_delta(w_rot, w_orb, obliq, i, sol_phase, orb_pos, sub_long, time, theta_1, phi)
array_2 = analytic_delta(w_rot, w_orb, obliq, i, sol_phase, orb_pos, sub_long, time, theta_2, phi)
array_3 = analytic_delta(w_rot, w_orb, obliq, i, sol_phase, orb_pos, sub_long, time, theta_3, phi)
array_4 = analytic_delta(w_rot, w_orb, obliq, i, sol_phase, orb_pos, sub_long, time, theta_4, phi)

plt.subplot(4,1,1)
plt.plot(array_1)
plt.ylabel('flux')
plt.xlabel('time(rotations)')

plt.subplot(4,1,2)
plt.plot(array_2)
plt.ylabel('flux')
plt.xlabel('time(rotations)')

plt.subplot(4,1,3)
plt.plot(array_3)
plt.ylabel('flux')
plt.xlabel('time(rotations)')


plt.subplot(4,1,4)
plt.plot(array_4)
plt.ylabel('flux')
plt.xlabel('time(rotations)')

plt.show()


