import numpy as np 

def analytic_delta(lat, long, w_rot, w_orb, obliq, i, sol_phase, orb_pos, sub_long, time, theta, phi): 

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
        print("All is well.")
    
    flux = []

    for t in time: 

        cos_phi_knot = np.cos(w_rot*t)
        #print("Debugging: cos_phi_knot:")
        #print(cos_phi_knot)

        cos_theta_knot = np.cos(i)*np.cos(obliq) + np.sin(i)*np.sin(obliq)*np.cos(sol_phase)
        #print("Debugging: cos_theta_knot:")
        #print(cos_theta_knot)

        sin_theta_knot = np.sqrt(1-(cos_theta_knot)**2)
        #print("Debugging: sin_theta_knot")
        #print(sin_theta_knot)

        sin_theta = np.sin(theta)
        #print("Debugging: sin_theta")
        #print(sin_theta)

        cos_theta = np.cos(theta)
        #print("Debugging: cos_theta")
        #print(cos_theta)

        cos_theta_s = np.sin(obliq)*np.cos(w_orb*t - sub_long)
        #print("Debugging: cos_theta_s")
        #print(cos_theta_s) 

        sin_phi_knot = -np.sin(w_rot*t)  
        #print("Debugging: sin_phi_knot")
        #print(sin_phi_knot) 

        sin_theta_s = np.sqrt(1-(cos_theta_s)**2)          
        #print("Debugging: sin_theta_s")
        #print(sin_theta_s)  

        cos_piw = np.sin(i)*np.cos(w_orb*t)
        #print("Debugging: cos_piw")
        #print(cos_piw) 

        cos_phi_knot_s = (cos_piw - (cos_theta_knot*cos_theta_s))/(sin_theta_knot * sin_theta_s)
        #print("Debugging: cos_phi_knot_s")
        #print(cos_phi_knot_s) 

        sin_phi_knot_s = np.sqrt(1-(cos_phi_knot_s)**2)                           
        #print("Debugging: sin_phi_knot_s:")
        #print(sin_phi_knot_s) 

        sin_phi_s = sin_phi_knot*cos_phi_knot_s + cos_phi_knot*sin_phi_knot_s
        #print("Debugging sin_phi_s:")
        #print(sin_phi_s)

        cos_phi_s = cos_phi_knot*cos_phi_knot_s - sin_phi_knot*sin_phi_knot_s
        #print("Debugging cos_phi_s")
        #print(cos_phi_s)

        #Visibility and illumination 

        cos_phi_phi_knot = np.cos(phi)*cos_phi_knot + np.sin(phi)*sin_phi_knot
        #print("Debugging cos_phi_phi_knot")
        #print(cos_phi_phi_knot)

        cos_phi_phi_s = np.cos(phi)*cos_phi_s + np.sin(phi)*sin_phi_s 
        #print("Debugging cos_phi_phi_s")
        #print(cos_phi_phi_s)

        Vnz = np.sin(theta)*sin_theta_knot*cos_phi_phi_knot + np.cos(theta)*cos_theta_knot 
        
        #print("Debugging Visbility nz: ")
        #print(Vnz)

        Inz = np.sin(theta)*sin_theta_s * cos_phi_phi_s + np.cos(theta)*cos_theta_s 
        #print("Debugging Illumination nz")
        #print(Inz)

        single_flux = (1/np.pi)*(max(Vnz, 0))*(max(Inz, 0))
        #print("Debugging single_flux")
        #print(single_flux)

        flux.append(single_flux)
    
    return flux 



#Testing
time = [1,1,1]
lat = 2
long = 1 
w_rot = 1
w_orb = 1
obliq = 1
i =1 
sol_phase = 1
orb_pos = 1
sub_long = 1
theta = np.pi/2
phi = np.pi/2

print(analytic_delta(lat, long, w_rot, w_orb, obliq, i, sol_phase, orb_pos, sub_long, time, theta, phi))


        




