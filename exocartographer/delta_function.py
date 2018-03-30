#STILL WORK IN PROGRESS

#importing modules 

import numpy as np
    

def analytic_delta(lat, long, w_rot, w_orb, obliq, i, solstice_phase, time):

    #First: checking to ensure that all parameters are in the proper values 


    #Ranges

    #w_orb 
    if (not(w_orb > 0)): 
        raise ValueError("Orbital angular frequency must be in ]0,+infinity[.")

    #obliq 
    elif (not((0 <= obliq) and ((np.pi*(1/2))>= obliq))): 
        raise ValueError("Obliquity angular frequency must be in [0,pi/2]")
    #i
    elif (not((i >= 0) and (i <= np.pi))): 
        raise ValueError("Orbital inclination must be between [0, pi]")

    #solstice phase 
    elif(not((solstice_phase >= 0) and (solstice_phase < (2*np.pi)))): 
        raise ValueError("Solstice phase must be in [0,2pi[")
    else: 
        print("All is well")

    #Setting each of the starting orbital positions for each value of t in the array 
    flux = []
    for i in time: 
        t = time[i]
        
        #print("Loop iteration ", t)
        starting_orbital_position = -1*w_rot*t
        sub_observer_longitude =  w_orb*t

        cos_theta_knot = np.cos(i)*np.cos(obliq) + np.sin(i)*np.sin(obliq)*np.cos(solstice_phase)
        sin_theta_knot = np.sqrt(1-(cos_theta_knot)**2)
        sin_theta_s = np.sqrt(1-np.sin(obliq)**2 * np.cos(w_orb*t - solstice_phase))
        cos_theta_s = np.sqrt(np.sin(obliq)*np.cos(w_orb*t - solstice_phase))

        #Debugging
        #print(cos_theta_knot)
        #print(sin_theta_knot)
        #print(sin_theta_s)
        #print(cos_theta_s)

        cos_piw = np.sin(i)*np.cos(w_orb*t)
        cos_sknot = ((cos_piw)-(cos_theta_knot*cos_theta_s))/(sin_theta_knot*sin_theta_s)

        sin_phi_s_phi_knot = (1/(sin_theta_knot*sin_theta_s))*(np.sqrt(1-(cos_theta_knot)**2 - (cos_theta_s)**2 - (cos_piw)**2 + 2*cos_theta_knot*cos_theta_s * cos_piw))
        cos_phi_s = np.cos(w_rot*t)*cos_sknot - ((-1)*np.sin(w_rot*t)*sin_phi_s_phi_knot)

        sin_phi_s = (((-np.sin(w_rot*t)*(np.sin(i)*np.cos(w_orb*t)- cos_theta_knot*np.sin(obliq)*np.cos(w_orb*t - solstice_phase)))
            + np.cos(w_rot*t)*(np.sin(i)*np.sin(w_orb*t)*np.cos(obliq)-np.cos(i)*np.sin(obliq)*np.sin(w_orb*t-solstice_phase)))/
            (np.sqrt(1-(np.cos(i)*np.cos(obliq) + np.sin(i)*np.cos(obliq)*np.cos(solstice_phase))**2) 
            * np.sqrt(1-((np.sin(obliq)**2)*(np.cos(w_orb*t - solstice_phase))**2))))

        cos_phi_phiknot = np.cos(lat)*np.cos(w_rot*t) + np.sin(lat)*((-1)*np.sin(w_rot*t))
        cos_phi_phis = np.cos(lat)*cos_phi_s + np.sin(lat)*sin_phi_s 

        #Illumination 
        I = np.sin(long)*sin_theta_s*cos_phi_phis + np.cos(long)*cos_theta_s
        #Visibility 
        v = np.sin(long)*sin_theta_knot*cos_phi_phiknot + np.cos(long)*cos_theta_knot

        flux.append((1/np.pi)*(max(v,0)*max((I,0)))) 

    return flux

#Testing
times = [1,2,3,4,5,5]
print(analytic_delta(1,1,1,1,1,1,1,times))


        




