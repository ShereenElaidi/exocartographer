import numpy as np, quaternion
import spherical_functions as sf


# a function to obtain the wigner D matrix for alpha,beta, gamma, ell, mp, and m
# values

def get_wigner_d(alpha, beta, gamma, ell, mp, m):
    R = quaternion.from_euler_angles(alpha, beta, gamma)
    wigner_d = sf.Wigner_D_element(R, ell, mp, m)
    return wigner_d

def get_phi(w, m):
    phi = 0
    if m ==0:
        phi = (1/2)*(np.sin(w) - w*np.cos(w))
    elif m == 2:
        phi = (1/4)*np.e**(1j*w)*(w-np.cos(w)*np.sin(w))
    elif m == -2:
        phi = (1/4)*np.e**(-1j*w)*(w-np.cos(w)*np.sin(w))
    else:
        phi = (2j*np.cos*w(1-np.e**(1j*m*w))-m*np.sin(w)*(1+np.e**(1j*m*w)))/(m*(m**2-4))
    return phi


def get_euler_angles(w, w_rot, w_orb, obliq, i, sol_phase, sub_long):

    # Making sure all input parameters are in bounds
    if(w_orb<0):
        raise ValueError("Values for w_orb must be bigger than 0.")
    elif (obliq > np.pi / 2 or obliq < 0):
        raise ValueError("Values for obliquity must be [0, pi/2]")
    # Inclination
    elif (i < 0 or i > np.pi):
        raise ValueError("Values for inclination must be [0,pi]")
    # Solstice phase
    elif (sol_phase < 0 or sol_phase >= 2 * np.pi):
        raise ValueError("Values for solstice phase must be [0,2pi[")

    cos_theta_knot = np.cos(i) * np.cos(obliq) + np.sin(i) * np.sin(obliq) * np.cos(sol_phase)
    sin_theta_s = np.sqrt(1-(cos_theta_s)**2)
    cos_phi_knot =  np.cos(w_rot*t)
    cos_theta_s = np.sin(obliq)*np.cos(w_orb*t - sub_long)
    sin_theta_knot = np.sqrt(1-(cos_theta_knot)**2)
    cos_phi_s = cos_phi_knot*cos_phi_knot_s - sin_phi_knot*sin_phi_knot_s
    sin_phi_knot  = -np.sin(w_rot*t)
    sin_phi_knot_s = np.sqrt(1-(cos_phi_knot_s)**2)
    sin_w = np.sin(w)
    cos_w = np.cos(w)
    sin_phi_s = sin_phi_knot*cos_phi_knot_s + cos_phi_knot*sin_phi_knot_s

    z_ell_2 = cos_theta_knot*sin_theta_s * cos_phi_s - cos_theta_s *sin_theta_knot*cos_theta_knot
    z_ell_1 = cos_theta_s*sin_theta_knot* sin_phi_knot - cos_theta_knot*sin_theta_s*sin_phi_s
    z_ell_3 = (sin_theta_knot * sin_theta_s*sin_phi_knot_s)/(sin_w)
    y_ell_3 = cos_theta_s*sin(w)
    x_ell_3 = cos_theta_knot+cos_theta_s*cos_w

    tan_alpha = z_ell_2/z_ell_1
    cos_beta = z_ell_3
    tan_gamma = y_ell_3/-x_ell_3

    alpha = math.atan2(tan_alpha)
    beta = np.arccos(cos_beta)
    gamma = math.atan2(tan_gamma)

    return alpha, beta, gamma
