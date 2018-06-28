from analytic_kernel import viewgeom
from analytic_kernel import kernel
import numpy as np
from scipy.integrate import dblquad

# TODO: ask about this:
def map(lat, lng, colat, colng, size):
    '''A function to mimic the delta function for the numerical solution. If the given longitude and latitude
    are within the bright pixel (defined by the size parameter), return 1. Else, return 0.

    Inputs:
    :param lat: the latitude of the point
    :param lng: the longitude of the point
    :param colat: the latitude of the bright point
    :param colng: the longitude of the bright point
    :param size: the size of the pixel.

    Output: 0 or 1 based on if the lat and lng fall within the pixel containing colat and colng'''

    if (lat <= size*colat or lat >= -1*size*colat) and (lng <= size*colng or lat >= -1*size*colng):
        return 1
    else:
        return 0

def get_kernel(times, wrot, worb, obq, inc, xisol, xi0, lat, lng):
    '''A function to obtain the kernel for the numerical integration solution.

    Inputs: an array of time-values (discrete), wrot, worb, obq, inc, xisol, xi0, lat, and lng

    Output: the value of the kernel'''
    view_geom = viewgeom(times, wrot, worb, obq, inc, xisol, xi0)
    sth = np.sin(lat)
    cth = np.cos(lat)
    sph = np.sin(lat)
    cph = np.cos(lng)
    return kernel(sth, cth, sph, cph, view_geom)

def numerically_integrate(lat, lng, times, wrot, worb, obq, inc, xisol, xi0, colat, colng, size):
    '''A function to implement the numerical solution: numerically integrating the varying flux from a small
    pixel. '''
    flux = []
    integrand = get_kernel(times, wrot, worb, obq, inc, xisol, xi0, lat, lng)*map(lat, lng, colat, colng, size)*np.sin(lat)
    for t in range(0, len(integrand)):
        flux.append(dblquad(lambda lat, lng: integrand[t], 0, np.pi, lambda theta:0,
                       lambda theta: 2*np.pi)[0])
    return flux