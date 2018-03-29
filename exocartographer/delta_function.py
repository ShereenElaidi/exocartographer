#STILL WORK IN PROGRESS

#importing modules 

import numpy as np

def analytic_delta(lat, long, w_rot, w_orb, obliq, i, solstice_phase, *time):

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

    for x in range(len(time)): 
        starting_orbital_position = -1*w_rot*time[x]
        sub_observer_longitude =  w_orb*time[x]

#Testing
analytic_delta(1,1,1,1,1,1,1,[1,1,1,1])


        




