import matplotlib
import numpy as np
import matplotlib.cm as cm
import matplotlib.mlab as mlab
import matplotlib.pyplot as plt
import pdb

# Obtaining Vnz and Inz
colat = np.array([np.pi, 3*(np.pi/4), np.pi/2, np.pi/4, 0])
lng = np.array([-np.pi, -np.pi/2, 0, np.pi/2, np.pi])


def Vnz (colat, lng, i, obq, sol_phase, times):
    times = times
    p_rotation = 23.934
    p_orbit = 365.256363 * 24.0
    w_rot = 2 * np.pi / p_rotation
    w_orb = 2 * np.pi / p_orbit
    i = i
    obq = obq
    sol_phase = sol_phase

    # ------------------
    theta = colat
    phi = lng

    # Appendix C: Derivation of the Time-Varying Geometry
    # (C1)
    cos_theta_knot = (np.cos(i) * np.cos(obq)) + (np.sin(i) * np.sin(obq) * np.cos(sol_phase))
    sin_theta_knot = np.sqrt(1 - (cos_theta_knot ** 2))

    # (C2)
    cos_phi_knot = np.cos(w_rot * times)
    sin_phi_knot = -np.sin(w_rot * times)

    # (C3)
    cos_theta_s = np.sin(obq) * np.cos((w_orb * times) - sol_phase)
    sin_theta_s = np.sqrt(1 - (((np.sin(obq)) ** 2) * ((np.cos((w_orb * times) - sol_phase)) ** 2)))
    cos_phi_s = 0
    sin_phi_s = 0

    if i == np.pi / 2:
        # (C1) Orbital Inclination
        cos_theta_knot = np.sin(obq) * np.cos(sol_phase)

        # (C12)
        cos_phi_s_1 = (np.cos(w_rot * times)) * (
                np.cos(w_orb * times) - cos_theta_knot * np.sin(obq) * np.cos(w_orb * times - sol_phase))
        cos_phi_s_2 = np.sin(w_rot * times) * np.sin(w_orb * times) * np.cos(obq)
        cos_phi_s_3 = sin_theta_knot * (
            np.sqrt(1 - ((((np.sin(obq)) ** 2) * ((np.cos(w_orb * times - sol_phase)) ** 2)))))
        cos_phi_s = (cos_phi_s_1 + cos_phi_s_2) / (cos_phi_s_3)

        # (C13)
        sin_phi_s_1 = -np.sin(w_rot * times) * (
                np.cos(w_orb * times) - cos_theta_knot * np.sin(obq) * np.cos(w_orb * times - sol_phase))
        sin_phi_s_2 = np.cos(w_rot * times) * np.sin(w_orb * times) * np.cos(obq)
        sin_phi_s_3 = sin_theta_knot * (
            np.sqrt(1 - (((np.sin(obq)) ** 2) * ((np.cos(w_orb * times - sol_phase)) ** 2))))
        sin_phi_s = (sin_phi_s_1 + sin_phi_s_2) / sin_phi_s_3

    elif i == 0:
        cos_theta_knot = np.cos(obq)

        # (C14)
        cos_phi_s_1 = -np.cos(w_rot * times) * np.cos(obq) * np.cos(w_orb * times - sol_phase)
        cos_phi_s_2 = np.sin(w_rot * times) * np.sin(w_orb * times - sol_phase)
        cos_phi_s_3 = np.sqrt(1 - ((((np.sin(obq)) ** 2)) * ((np.cos(w_orb * times - sol_phase)) ** 2)))
        cos_phi_s = (cos_phi_s_1 - cos_phi_s_2) / cos_phi_s_3

        # (C15)
        sin_phi_s_1 = np.sin(w_rot * times) * np.cos(obq) * np.cos(w_orb * times - sol_phase)
        sin_phi_s_2 = np.cos(w_rot * times) * np.sin(w_orb * times - sol_phase)
        sin_phi_s_3 = np.sqrt(1 - (((np.sin(obq)) ** 2) * ((np.cos(w_orb * times - sol_phase)) ** 2)))
        sin_phi_s = (sin_phi_s_1 - sin_phi_s_2) / sin_phi_s_3

    elif obq == 0:

        # (C3) Planetary Obliquity
        theta_knot = i
        theta_s = np.pi / 2
        phi_s = (w_orb - w_rot) * times
        cos_theta_knot = np.cos(i)
        cos_theta_s = 0
        cos_phi_s = np.cos((w_orb - w_rot) * times)
        sin_phi_s = np.sin((w_orb - w_rot) * times)

    elif obq == np.pi / 2:
        # (C20)
        cos_theta_knot = np.sin(i) * np.sin(obq) * np.cos(sol_phase)
        sin_theta_knot = np.sqrt(1 - (cos_theta_knot ** 2))

        # (C21)
        cos_theta_s = np.cos(w_orb * times - sol_phase)
        sin_theta_s = np.sin(w_orb * times - sol_phase)

        # (C22)
        cos_phi_s_1 = np.cos(w_rot * times) * (
                np.sin(i) * np.cos(w_orb * times) - cos_theta_knot * np.cos(w_orb * times - sol_phase))
        cos_phi_s_2 = np.sin(w_rot * times) * np.cos(i) * np.sin(w_orb * times - sol_phase)
        cos_phi_s_3 = np.sqrt(1 - (((np.sin(i)) ** 2) * ((np.cos(sol_phase)) ** 2)))
        cos_phi_s_4 = np.sin(w_orb * times - sol_phase)
        cos_phi_s = (cos_phi_s_1 - cos_phi_s_2) / (cos_phi_s_3 * cos_phi_s_4)

        # (C23)
        sin_phi_s_1 = -np.sin(w_rot * times) * (
                np.sin(i) * np.cos(w_orb * times) - cos_theta_knot * np.cos(w_orb * times - sol_phase))
        sin_phi_s_2 = np.cos(w_rot * times) * np.cos(i) * np.sin(w_orb * times - sol_phase)
        sin_phi_s_3 = np.sqrt(1 - (((np.sin(i)) ** 2) * ((np.cos(sol_phase)) ** 2)))
        sin_phi_s_4 = np.sin(w_orb * times - sol_phase)
        sin_phi_s = (sin_phi_s_1 - sin_phi_s_2) / (sin_phi_s_3 * sin_phi_s_4)

    else:
        # (C4)
        cos_phi_s_1 = np.cos(w_rot * times) * (
                np.sin(i) * np.cos(w_orb * times) - cos_theta_knot * np.sin(obq) * np.cos(w_orb * times - sol_phase))
        cos_phi_s_2 = np.sin(w_rot * times) * (
                np.sin(i) * np.sin(w_orb * times) * np.cos(obq) - np.cos(i) * np.sin(obq) * np.sin(
            w_orb * times - sol_phase))
        cos_phi_s_3 = np.sqrt(1 - (np.cos(i) * np.cos(obq) + np.sin(i) * np.sin(obq) * np.cos(sol_phase)) ** 2)
        cos_phi_s_4 = np.sqrt(1 - (((np.sin(obq)) ** 2) * ((np.cos(w_orb * times - sol_phase)) ** 2)))
        cos_phi_s = (cos_phi_s_1 + cos_phi_s_2) / (cos_phi_s_3 * cos_phi_s_4)

        # (C5)
        sin_phi_s_1 = -np.sin(w_rot * times) * (
                np.sin(i) * np.cos(w_orb * times) - cos_theta_knot * np.sin(obq) * np.cos(w_orb * times - sol_phase))
        sin_phi_s_2 = np.cos(w_rot * times) * (
                np.sin(i) * np.sin(w_orb * times) * np.cos(obq) - np.cos(i) * np.sin(obq) * np.sin(
            w_orb * times - sol_phase))
        sin_phi_s_3 = np.sqrt(1 - (np.cos(i) * np.cos(obq) + np.sin(i) * np.sin(obq) * np.cos(sol_phase)) ** 2)
        sin_phi_s_4 = np.sqrt(1 - (((np.sin(obq)) ** 2) * ((np.cos(w_orb * times - sol_phase)) ** 2)))
        sin_phi_s = (sin_phi_s_1 + sin_phi_s_2) / (sin_phi_s_3 * sin_phi_s_4)

    # Delta Maps
    cos_phi_phi_s = np.cos(phi) * cos_phi_s + np.sin(phi) * sin_phi_s
    cos_phi_phi_knot = np.cos(phi) * cos_phi_knot + np.sin(phi) * sin_phi_knot

    # (4)
    Inz = np.sin(theta) * sin_theta_s * cos_phi_phi_s + np.cos(theta) * cos_theta_s

    # (5)
    Vnz = np.sin(theta) * sin_theta_knot * cos_phi_phi_knot + np.cos(theta) * cos_theta_knot
    return Vnz

def Inz (colat, lng, i, obq, sol_phase, times):
    times = times
    p_rotation = 23.934
    p_orbit = 365.256363 * 24.0
    w_rot = 2 * np.pi / p_rotation
    w_orb = 2 * np.pi / p_orbit
    i = i
    obq = obq
    sol_phase = sol_phase

    # ------------------
    theta = colat
    phi = lng

    # Appendix C: Derivation of the Time-Varying Geometry
    # (C1)
    cos_theta_knot = (np.cos(i) * np.cos(obq)) + (np.sin(i) * np.sin(obq) * np.cos(sol_phase))
    sin_theta_knot = np.sqrt(1 - (cos_theta_knot ** 2))

    # (C2)
    cos_phi_knot = np.cos(w_rot * times)
    sin_phi_knot = -np.sin(w_rot * times)

    # (C3)
    cos_theta_s = np.sin(obq) * np.cos((w_orb * times) - sol_phase)
    sin_theta_s = np.sqrt(1 - (((np.sin(obq)) ** 2) * ((np.cos((w_orb * times) - sol_phase)) ** 2)))
    cos_phi_s = 0
    sin_phi_s = 0

    if i == np.pi / 2:
        # (C1) Orbital Inclination
        cos_theta_knot = np.sin(obq) * np.cos(sol_phase)

        # (C12)
        cos_phi_s_1 = (np.cos(w_rot * times)) * (
                np.cos(w_orb * times) - cos_theta_knot * np.sin(obq) * np.cos(w_orb * times - sol_phase))
        cos_phi_s_2 = np.sin(w_rot * times) * np.sin(w_orb * times) * np.cos(obq)
        cos_phi_s_3 = sin_theta_knot * (
            np.sqrt(1 - ((((np.sin(obq)) ** 2) * ((np.cos(w_orb * times - sol_phase)) ** 2)))))
        cos_phi_s = (cos_phi_s_1 + cos_phi_s_2) / (cos_phi_s_3)

        # (C13)
        sin_phi_s_1 = -np.sin(w_rot * times) * (
                np.cos(w_orb * times) - cos_theta_knot * np.sin(obq) * np.cos(w_orb * times - sol_phase))
        sin_phi_s_2 = np.cos(w_rot * times) * np.sin(w_orb * times) * np.cos(obq)
        sin_phi_s_3 = sin_theta_knot * (
            np.sqrt(1 - (((np.sin(obq)) ** 2) * ((np.cos(w_orb * times - sol_phase)) ** 2))))
        sin_phi_s = (sin_phi_s_1 + sin_phi_s_2) / sin_phi_s_3

    elif i == 0:
        cos_theta_knot = np.cos(obq)

        # (C14)
        cos_phi_s_1 = -np.cos(w_rot * times) * np.cos(obq) * np.cos(w_orb * times - sol_phase)
        cos_phi_s_2 = np.sin(w_rot * times) * np.sin(w_orb * times - sol_phase)
        cos_phi_s_3 = np.sqrt(1 - ((((np.sin(obq)) ** 2)) * ((np.cos(w_orb * times - sol_phase)) ** 2)))
        cos_phi_s = (cos_phi_s_1 - cos_phi_s_2) / cos_phi_s_3

        # (C15)
        sin_phi_s_1 = np.sin(w_rot * times) * np.cos(obq) * np.cos(w_orb * times - sol_phase)
        sin_phi_s_2 = np.cos(w_rot * times) * np.sin(w_orb * times - sol_phase)
        sin_phi_s_3 = np.sqrt(1 - (((np.sin(obq)) ** 2) * ((np.cos(w_orb * times - sol_phase)) ** 2)))
        sin_phi_s = (sin_phi_s_1 - sin_phi_s_2) / sin_phi_s_3

    elif obq == 0:

        # (C3) Planetary Obliquity
        theta_knot = i
        theta_s = np.pi / 2
        phi_s = (w_orb - w_rot) * times
        cos_theta_knot = np.cos(i)
        cos_theta_s = 0
        cos_phi_s = np.cos((w_orb - w_rot) * times)
        sin_phi_s = np.sin((w_orb - w_rot) * times)

    elif obq == np.pi / 2:
        # (C20)
        cos_theta_knot = np.sin(i) * np.sin(obq) * np.cos(sol_phase)
        sin_theta_knot = np.sqrt(1 - (cos_theta_knot ** 2))

        # (C21)
        cos_theta_s = np.cos(w_orb * times - sol_phase)
        sin_theta_s = np.sin(w_orb * times - sol_phase)

        # (C22)
        cos_phi_s_1 = np.cos(w_rot * times) * (
                np.sin(i) * np.cos(w_orb * times) - cos_theta_knot * np.cos(w_orb * times - sol_phase))
        cos_phi_s_2 = np.sin(w_rot * times) * np.cos(i) * np.sin(w_orb * times - sol_phase)
        cos_phi_s_3 = np.sqrt(1 - (((np.sin(i)) ** 2) * ((np.cos(sol_phase)) ** 2)))
        cos_phi_s_4 = np.sin(w_orb * times - sol_phase)
        cos_phi_s = (cos_phi_s_1 - cos_phi_s_2) / (cos_phi_s_3 * cos_phi_s_4)

        # (C23)
        sin_phi_s_1 = -np.sin(w_rot * times) * (
                np.sin(i) * np.cos(w_orb * times) - cos_theta_knot * np.cos(w_orb * times - sol_phase))
        sin_phi_s_2 = np.cos(w_rot * times) * np.cos(i) * np.sin(w_orb * times - sol_phase)
        sin_phi_s_3 = np.sqrt(1 - (((np.sin(i)) ** 2) * ((np.cos(sol_phase)) ** 2)))
        sin_phi_s_4 = np.sin(w_orb * times - sol_phase)
        sin_phi_s = (sin_phi_s_1 - sin_phi_s_2) / (sin_phi_s_3 * sin_phi_s_4)

    else:
        # (C4)
        cos_phi_s_1 = np.cos(w_rot * times) * (
                np.sin(i) * np.cos(w_orb * times) - cos_theta_knot * np.sin(obq) * np.cos(w_orb * times - sol_phase))
        cos_phi_s_2 = np.sin(w_rot * times) * (
                np.sin(i) * np.sin(w_orb * times) * np.cos(obq) - np.cos(i) * np.sin(obq) * np.sin(
            w_orb * times - sol_phase))
        cos_phi_s_3 = np.sqrt(1 - (np.cos(i) * np.cos(obq) + np.sin(i) * np.sin(obq) * np.cos(sol_phase)) ** 2)
        cos_phi_s_4 = np.sqrt(1 - (((np.sin(obq)) ** 2) * ((np.cos(w_orb * times - sol_phase)) ** 2)))
        cos_phi_s = (cos_phi_s_1 + cos_phi_s_2) / (cos_phi_s_3 * cos_phi_s_4)

        # (C5)
        sin_phi_s_1 = -np.sin(w_rot * times) * (
                np.sin(i) * np.cos(w_orb * times) - cos_theta_knot * np.sin(obq) * np.cos(w_orb * times - sol_phase))
        sin_phi_s_2 = np.cos(w_rot * times) * (
                np.sin(i) * np.sin(w_orb * times) * np.cos(obq) - np.cos(i) * np.sin(obq) * np.sin(
            w_orb * times - sol_phase))
        sin_phi_s_3 = np.sqrt(1 - (np.cos(i) * np.cos(obq) + np.sin(i) * np.sin(obq) * np.cos(sol_phase)) ** 2)
        sin_phi_s_4 = np.sqrt(1 - (((np.sin(obq)) ** 2) * ((np.cos(w_orb * times - sol_phase)) ** 2)))
        sin_phi_s = (sin_phi_s_1 + sin_phi_s_2) / (sin_phi_s_3 * sin_phi_s_4)

    # Delta Maps
    cos_phi_phi_s = np.cos(phi) * cos_phi_s + np.sin(phi) * sin_phi_s
    cos_phi_phi_knot = np.cos(phi) * cos_phi_knot + np.sin(phi) * sin_phi_knot

    # (4)
    Inz = np.sin(theta) * sin_theta_s * cos_phi_phi_s + np.cos(theta) * cos_theta_s

    # (5)
    Vnz = np.sin(theta) * sin_theta_knot * cos_phi_phi_knot + np.cos(theta) * cos_theta_knot
    return Inz



#----------------------------------------------

X,Y = np.meshgrid(lng, colat)
colat = X
lng = Y
i =0
obq = 0
sol_phase = 0
times=0
f1 = Inz(colat=colat, lng=lng, i=i, obq=obq, sol_phase=sol_phase, times=times)
f2 = Vnz(colat=colat, lng=lng, i=i, obq=obq, sol_phase=sol_phase, times=times)
print("Illumination Nz")
print(f1)
print("Visibility Nz")
print(f2)
f3 = np.multiply(f1, f2)
fig = plt.figure()
ax = fig.add_subplot(111)
cpf = ax.contourf(X,Y,f3, 20, cmap=cm.Greys_r)
print(cpf.levels)
colours = ['k' if level<0 else 'w' for level in cpf.levels]
cp = ax.contour(X,Y,f3, colors=colours)
ax.clabel(cp, fontsize=12, colors=colours)
plt.contour(X,Y,f1, colors='red')
plt.contour(X,Y,f2, colors='blue')
plt.contour(X,Y,f3, colors='green')
plt.ylim(np.pi, 0)
plt.xlim(-np.pi, np.pi)
plt.show()

