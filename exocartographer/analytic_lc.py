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

    print("Inclination pi/2?")
    print(i==np.pi/2)
    print("Inclination 0?")
    print(i==0)
    print("Obliquity 0?")
    print(obq == 0)
    print("Obliquity pi/2?")
    print(obq == np.pi/2)

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
    cos_phi_phi_s = np.cos(phi)*cos_phi_s + np.sin(phi)*sin_phi_s
    cos_phi_phi_knot = np.cos(phi)*cos_phi_knot + np.sin(phi)*sin_phi_knot

    # (4)
    Inz = np.sin(theta)*sin_theta_s*cos_phi_phi_s + np.cos(theta)*cos_theta_s

    # (5)
    Vnz = np.sin(theta)*sin_theta_knot*cos_phi_phi_knot + np.cos(theta)*cos_theta_knot

    # (2)
    z = np.zeros_like(Inz)
    kernel = (1/np.pi)*(np.maximum(Vnz, z))*(np.maximum(Inz, z))
    # degrees_to_str = (1/3282.8)
    domega = hp.nside2pixarea(nside=nside, degrees=False)*np.sin(theta)
    return kernel*domega