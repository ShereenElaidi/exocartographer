from analytic_kernel import viewgeom
from analytic_kernel import kernel
import numpy as np
from scipy.integrate import dblquad


def map(lat, lng):
    return lat*lng

def get_kernel(lat, lng, t, wrot, worb, obq, inc, xisol, xi0):
    view = viewgeom(t,wrot, worb, obq, inc, xisol, xi0)
    sth = np.sin(lat)
    cth = np.cos(lat)
    sph = np.sin(lng)
    cph = np.cos(lng)
    return kernel(sth, cth,sph, cph, view)


lat = 0
lng = 0
worb = 2
t = np.array([1])
wrot = 1
orb = 1
obq = np.pi/3
inc = np.pi/2
xisol = np.pi
xi0 = 1

if (worb < 0):
    raise ValueError("Values for w_orb must be bigger than 0.")
# Obliquity
elif (obq > np.pi / 2 or obq < 0):
    raise ValueError("Values for obliquity must be [0, pi/2]")
# Inclination
elif (inc < 0 or inc > np.pi):
    raise ValueError("Values for inclination must be [0,pi]")
# Solstice phase
elif (xisol < 0 or xisol >= 2 * np.pi):
    raise ValueError("Values for solstice phase must be [0,2pi[")
else:
    print("")

def compute_flux(lat,lng, t, wrot, worb, obq, inc, xisol, xi0):
    integrand = get_kernel(lat, lng, t, wrot, worb, obq, inc, xisol, xi0) * map(lat, lng) * np.sin(lat)
    flux = dblquad(lambda lat, lng: integrand, 0, np.pi, lambda theta: 0, lambda theta: 2 * np.pi)[0]
    print(flux)

