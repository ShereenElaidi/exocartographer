import numpy as np

def analytic_delta(w_rot, w_orb, obliq, i, sol_phase, orb_pos, sub_long, time, theta, phi): 

    
    flux = []

    for t in time: 
        cos_phi_knot = np.cos(w_rot*t)
        cos_theta_knot = np.cos(i)*np.cos(obliq) + np.sin(i)*np.sin(obliq)*np.cos(sol_phase)
        sin_theta_knot = np.sqrt(1-(cos_theta_knot)**2)
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