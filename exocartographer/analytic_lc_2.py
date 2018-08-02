import numpy as np
import pdb


def analytic_lc_2(pt_lat, pt_lng, lat, lng, w_rot, w_orb, obq, i, sol_phase, times, nside):

    if (pt_lat != lat) or (pt_lng != lng):
        return np.zeros_like(times)

    domega = 1

    theta = lat
    phi = lng

    cos_theta_knot = (np.cos(i)*np.cos(obq)) + (np.sin(i)*np.sin(obq)*np.cos(sol_phase))
    sin_theta_knot = np.sqrt(1-(cos_theta_knot**2))

    cos_phi_knot = np.cos(w_rot*times)
    sin_phi_knot = -np.sin(w_rot*times)

    cos_theta_s = np.sin(obq)*np.cos(w_orb*times-sol_phase)
    sin_theta_s = np.sqrt(1-(((np.sin(obq))**2)*((np.cos(w_orb*times-sol_phase))**2)))


    # cos_phi_s
    cos_phi_s_1 = (np.cos(w_rot*times))*(np.sin(i)*np.cos(w_orb*times) - cos_theta_knot*np.sin(obq)*np.cos(w_orb*times-sol_phase))
    cos_phi_s_2 = (np.sin(w_rot*times))*(np.sin(i)*np.sin(w_orb*times)*np.cos(obq) - np.cos(i)*np.sin(obq)*np.sin(w_orb*times-sol_phase))
    cos_phi_s_3 = np.sqrt(1-(np.cos(i)*np.cos(obq) + np.sin(i)*np.sin(obq)*np.cos(sol_phase))**2)
    cos_phi_s_4 = np.sqrt(1-(((np.sin(obq))**2)*((np.cos(w_orb*times-sol_phase))**2)))
    cos_phi_s = (cos_phi_s_1+cos_phi_s_2)/(cos_phi_s_3*cos_phi_s_4)

    # sin_phi_s
    sin_phi_s_1 = (-np.sin(w_rot*times))*(np.sin(i)*np.cos(w_orb*times) - cos_theta_knot*np.sin(obq)*np.cos(w_orb*times-sol_phase))
    sin_phi_s_2 = (np.cos(w_rot*times))*(np.sin(i)*np.sin(w_orb*times)*np.cos(obq) - np.cos(i)*np.sin(obq)*np.sin(w_orb*times-sol_phase))
    sin_phi_s_3 = cos_phi_s_3
    sin_phi_s_4 = cos_phi_s_4
    sin_phi_s = (sin_phi_s_1 + sin_phi_s_2)/(sin_phi_s_3*sin_phi_s_4)


    # cos_phi_phi_knot
    cos_phi_phi_knot = np.cos(phi)*cos_phi_knot + np.sin(phi)*sin_phi_knot

    # cos_phi_phi_s
    cos_phi_phi_s = np.cos(phi)*cos_phi_s + np.sin(phi)*sin_phi_s

    Inz = np.sin(theta)*sin_theta_s*cos_phi_phi_s + np.cos(theta)*cos_theta_s
    Vnz = np.sin(theta)*sin_theta_knot*cos_phi_phi_knot+ np.cos(theta)*cos_theta_knot
    z = np.zeros(Inz.size)
    return domega * (1 / np.pi) * (np.maximum(Vnz, z)) * (np.maximum(Inz, z))






