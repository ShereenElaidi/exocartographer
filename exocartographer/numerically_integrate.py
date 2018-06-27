from analytic_kernel import viewgeom
from analytic_kernel import kernel
import numpy as np
from scipy.integrate import dblquad

# Basic integration: double integral, Kernel*map*sintheta


# TODO: ask about this:
def map(lat, lng):
    return lat*lng

def get_kernel(times, wrot, worb, obq, inc, xisol, xi0, lat, lng):
    view_geom = viewgeom(times, wrot, worb, obq, inc, xisol, xi0)
    sth = np.sin(lat)
    cth = np.cos(lat)
    sph = np.sin(lat)
    cph = np.cos(lng)
    return kernel(sth, cth, sph, cph, view_geom)

def numerically_integrate(lat, lng, times, wrot, worb, obq, inc, xisol, xi0):
    flux = []
    integrand = get_kernel(times, wrot, worb, obq, inc, xisol, xi0, lat, lng)
    for t in range(0, len(integrand)):
        flux.append(dblquad(lambda lat, lng: integrand[t], 0, np.pi, lambda theta:0,
                       lambda theta: 2*np.pi)[0])
    return flux