import numpy as np 

def analytic_delta(w_rot, w_orb, obliq, i, sol_phase, orb_pos, sub_long, time, theta, phi): 

    #First checking to make sure all input parameters are in bounds 

    #w_orb
    if(w_orb < 0):
        raise ValueError("Values for w_orb must be bigger than 0.")

    #Obliquity 
    elif (obliq > np.pi/2 or obliq < 0 ): 
        raise ValueError("Values for obliquity must be [0, pi/2]")

    #Inclination 
    elif (i < 0 or i > np.pi): 
        raise ValueError("Values for inclination must be [0,pi]")
    
    #Solstice phase 
    elif (sol_phase < 0 or sol_phase >= 2*np.pi): 
        raise ValueError("Values for solstice phase must be [0,2pi[")
    else: 
        print("")
    
    flux = []

    for t in time: 
        cos_phi_knot = np.cos(w_rot*t)
        cos_theta_knot = np.cos(i)*np.cos(obliq) + np.sin(i)*np.sin(obliq)*np.cos(sol_phase)
        sin_theta_knot = np.sqrt(1-(cos_theta_knot)**2)
        sin_theta = np.sin(theta)
        cos_theta = np.cos(theta)
        cos_theta_s = np.sin(obliq)*np.cos(w_orb*t - sub_long)
        sin_phi_knot = -np.sin(w_rot*t)  
        sin_theta_s = np.sqrt(1-(cos_theta_s)**2)          
        cos_piw = np.sin(i)*np.cos(w_orb*t)
        cos_phi_knot_s = (cos_piw - (cos_theta_knot*cos_theta_s))/(sin_theta_knot * sin_theta_s)
        sin_phi_knot_s = np.sqrt(1-(cos_phi_knot_s)**2)                           
        sin_phi_s = sin_phi_knot*cos_phi_knot_s + cos_phi_knot*sin_phi_knot_s
        cos_phi_s = cos_phi_knot*cos_phi_knot_s - sin_phi_knot*sin_phi_knot_s

        #Visibility and illumination 
        cos_phi_phi_knot = np.cos(phi)*cos_phi_knot + np.sin(phi)*sin_phi_knot
        cos_phi_phi_s = np.cos(phi)*cos_phi_s + np.sin(phi)*sin_phi_s 
        Vnz = np.sin(theta)*sin_theta_knot*cos_phi_phi_knot + np.cos(theta)*cos_theta_knot 
        Inz = np.sin(theta)*sin_theta_s * cos_phi_phi_s + np.cos(theta)*cos_theta_s 
        single_flux = (1/np.pi)*(max(Vnz, 0))*(max(Inz, 0))
        flux.append(single_flux)
    return flux 

#Testing
time = []
y = 0
for y in range (0, 100):
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


import matplotlib.pyplot as plt
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
exit(0)
        




