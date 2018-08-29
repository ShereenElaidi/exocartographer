import numpy as np
import healpy as hp
import pdb

def analytic_lightcurve_2(colat, lng, w_rot, w_orb, obq, i, sol_phase, times, nside):
    # Function takes as input the co-latitude and the longitude
    theta = colat
    phi = lng

    # Appendix C: Derivation of the Time-Varying Geometry
    # (C1)
    cos_theta_knot = (np.cos(i)*np.cos(obq)) + (np.sin(i)*np.sin(obq)*np.cos(sol_phase))
    sin_theta_knot = np.sqrt(1-(cos_theta_knot**2))

    # (C2)
    cos_phi_knot = np.cos(w_rot*times)
    sin_phi_knot = -np.sin(w_rot*times)

    # (C3)
    cos_theta_s = np.sin(obq)*np.cos((w_orb*times)-sol_phase)
    sin_theta_s = np.sqrt(1-(((np.sin(obq))**2)*((np.cos((w_orb*times)-sol_phase))**2)))
    cos_phi_s = 0
    sin_phi_s = 0

    if i == np.pi/2:
        # (C1) Orbital Inclination
        cos_theta_knot = np.sin(obq)*np.cos(sol_phase)

        # (C12)
        cos_phi_s_1 = (np.cos(w_rot*times))*(np.cos(w_orb*times) - cos_theta_knot*np.sin(obq)*np.cos(w_orb*times-sol_phase))
        cos_phi_s_2 = np.sin(w_rot*times)*np.sin(w_orb*times)*np.cos(obq)
        cos_phi_s_3 = sin_theta_knot*(np.sqrt(1-((((np.sin(obq))**2)*((np.cos(w_orb*times-sol_phase))**2)))))
        cos_phi_s = (cos_phi_s_1+cos_phi_s_2)/(cos_phi_s_3)

        # (C13)
        sin_phi_s_1 = -np.sin(w_rot*times)*(np.cos(w_orb*times)-cos_theta_knot*np.sin(obq)*np.cos(w_orb*times-sol_phase))
        sin_phi_s_2 = np.cos(w_rot*times)*np.sin(w_orb*times)*np.cos(obq)
        sin_phi_s_3 = sin_theta_knot*(np.sqrt(1-(((np.sin(obq))**2)*((np.cos(w_orb*times-sol_phase))**2))))
        sin_phi_s = (sin_phi_s_1+sin_phi_s_2)/sin_phi_s_3

    elif i == 0 :
        cos_theta_knot = np.cos(obq)

        # (C14)
        cos_phi_s_1 = -np.cos(w_rot*times)*np.cos(obq)*np.cos(w_orb*times-sol_phase)
        cos_phi_s_2 = np.sin(w_rot*times)*np.sin(w_orb*times-sol_phase)
        cos_phi_s_3 = np.sqrt(1-((((np.sin(obq))**2))*((np.cos(w_orb*times-sol_phase))**2)))
        cos_phi_s = (cos_phi_s_1-cos_phi_s_2)/cos_phi_s_3

        # (C15)
        sin_phi_s_1 = np.sin(w_rot*times)*np.cos(obq)*np.cos(w_orb*times-sol_phase)
        sin_phi_s_2 = np.cos(w_rot*times)*np.sin(w_orb*times-sol_phase)
        sin_phi_s_3 = np.sqrt(1-(((np.sin(obq))**2)*((np.cos(w_orb*times-sol_phase))**2)))
        sin_phi_s = (sin_phi_s_1-sin_phi_s_2)/sin_phi_s_3

    elif obq == 0:

        # (C3) Planetary Obliquity
        theta_knot = i
        theta_s = np.pi/2
        phi_s = (w_orb-w_rot)*times
        cos_theta_knot = np.cos(i)
        cos_theta_s = 0
        cos_phi_s = np.cos((w_orb-w_rot)*times)
        sin_phi_s = np.sin((w_orb-w_rot)*times)

    elif obq == np.pi/2:
        # (C20)
        cos_theta_knot = np.sin(i)*np.sin(obq)*np.cos(sol_phase)
        sin_theta_knot = np.sqrt(1-(cos_theta_knot**2))

        # (C21)
        cos_theta_s = np.cos(w_orb*times-sol_phase)
        sin_theta_s = np.sin(w_orb*times-sol_phase)

        # (C22)
        cos_phi_s_1 = np.cos(w_rot*times)*(np.sin(i)*np.cos(w_orb*times) - cos_theta_knot*np.cos(w_orb*times-sol_phase))
        cos_phi_s_2 = np.sin(w_rot*times)*np.cos(i)*np.sin(w_orb*times-sol_phase)
        cos_phi_s_3 = np.sqrt(1-(((np.sin(i))**2)*((np.cos(sol_phase))**2)))
        cos_phi_s_4 = np.sin(w_orb*times-sol_phase)
        cos_phi_s = (cos_phi_s_1-cos_phi_s_2)/(cos_phi_s_3*cos_phi_s_4)

        # (C23)
        sin_phi_s_1 = -np.sin(w_rot*times)*(np.sin(i)*np.cos(w_orb*times) - cos_theta_knot*np.cos(w_orb*times-sol_phase))
        sin_phi_s_2 = np.cos(w_rot*times)*np.cos(i)*np.sin(w_orb*times-sol_phase)
        sin_phi_s_3 = np.sqrt(1-(((np.sin(i))**2)*((np.cos(sol_phase))**2)))
        sin_phi_s_4 = np.sin(w_orb*times-sol_phase)
        sin_phi_s = (sin_phi_s_1-sin_phi_s_2)/(sin_phi_s_3*sin_phi_s_4)

    else:
        # (C4)
        cos_phi_s_1 = np.cos(w_rot*times)*(np.sin(i)*np.cos(w_orb*times) - cos_theta_knot*np.sin(obq)*np.cos(w_orb*times-sol_phase))
        cos_phi_s_2 = np.sin(w_rot*times)*(np.sin(i)*np.sin(w_orb*times)*np.cos(obq) - np.cos(i)*np.sin(obq)*np.sin(w_orb*times-sol_phase))
        cos_phi_s_3 = np.sqrt(1-(np.cos(i)*np.cos(obq) + np.sin(i)*np.sin(obq)*np.cos(sol_phase))**2)
        cos_phi_s_4 = np.sqrt(1-(((np.sin(obq))**2)*((np.cos(w_orb*times-sol_phase))**2)))
        cos_phi_s = (cos_phi_s_1+cos_phi_s_2)/(cos_phi_s_3*cos_phi_s_4)

        # (C5)
        sin_phi_s_1 = -np.sin(w_rot*times)*(np.sin(i)*np.cos(w_orb*times) - cos_theta_knot*np.sin(obq)*np.cos(w_orb*times-sol_phase))
        sin_phi_s_2 = np.cos(w_rot*times)*(np.sin(i)*np.sin(w_orb*times)*np.cos(obq) - np.cos(i)*np.sin(obq)*np.sin(w_orb*times-sol_phase))
        sin_phi_s_3 = np.sqrt(1-(np.cos(i)*np.cos(obq)+np.sin(i)*np.sin(obq)*np.cos(sol_phase))**2)
        sin_phi_s_4 = np.sqrt(1-(((np.sin(obq))**2)*((np.cos(w_orb*times-sol_phase))**2)))
        sin_phi_s = (sin_phi_s_1+sin_phi_s_2)/(sin_phi_s_3*sin_phi_s_4)

    # Delta Maps
    s_phi_o = -np.sin(w_rot*times)
    phi_o = np.arcsin(s_phi_o)
    phi_s = np.arcsin(sin_phi_s)
    s_th_s = np.sqrt(1 - (((np.sin(obq)) ** 2) * ((np.cos((w_orb * times) - sol_phase)) ** 2)))
    th_s = np.arcsin(s_th_s)
    c_th_o = (np.cos(i)*np.cos(obq)) + (np.sin(i)*np.sin(obq)*np.cos(sol_phase))
    th_o = np.full_like(times, np.arccos(c_th_o))

    arr_theta = np.full_like(times, theta)
    arr_phi = np.full_like(times, phi)

    Vnz = (np.sin(theta) * np.sin(th_o) * np.cos(arr_phi - phi_o)) + (np.cos(arr_theta) * np.cos(th_o))
    Inz = (np.sin(theta) * np.sin(th_s) * np.cos(arr_phi - phi_s)) + (np.cos(arr_theta) * np.cos(th_s))

    # (2)
    Vnz[Vnz < 0] = 0
    Inz[Inz < 0] = 0
    return ((Vnz*Inz)/np.pi)