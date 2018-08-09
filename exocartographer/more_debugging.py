# Import statements
import numpy as np
import matplotlib.pyplot as plt
from astropy_healpix import HEALPix
from astropy.coordinates import Angle
import astropy.units as u
import matplotlib.gridspec  as gridspec
import exocartographer
import math
import healpy as hp

# First, some helper functions. cos_insolation is the equivalent to
# the illumination for the analytic solution. cos_obs is the equivalent
# to the visibility for the analytic solution.

p_rotation = 23.934
p_orbit = 365.256363 * 24.0
phi_orb = 0
inclination = np.pi/6

solstice_phase =  np.pi/6
obl_orientation = solstice_phase

obliquity = np.pi/6
phi_rot = solstice_phase

times=12

def _observer_normal_orbit_coords(inclination=inclination):
    cos_inc = np.cos(inclination)
    sin_inc = np.sqrt(1.0-cos_inc*cos_inc)
    return np.array([-sin_inc, 0.0, cos_inc])


def rotation_quaternions(axis, angles):
    angles = np.atleast_1d(angles)
    result = np.zeros((angles.shape[0], 4))
    result[:, 0] = np.cos(angles/2.0)
    result[:, 1:] = np.sin(angles/2.0)[:, np.newaxis]*axis
    return result

def quaternion_multiply(qa, qb):
    result = np.zeros(np.broadcast(qa, qb).shape)

    result[..., 0] = qa[..., 0]*qb[..., 0] - np.sum(qa[..., 1:]*qb[..., 1:], axis=-1)
    result[..., 1] = qa[..., 0]*qb[..., 1] + qa[..., 1]*qb[..., 0] + qa[..., 2]*qb[..., 3] - qa[..., 3]*qb[..., 2]
    result[..., 2] = qa[..., 0]*qb[..., 2] - qa[..., 1]*qb[..., 3] + qa[..., 2]*qb[..., 0] + qa[..., 3]*qb[..., 1]
    result[..., 3] = qa[..., 0]*qb[..., 3] + qa[..., 1]*qb[..., 2] - qa[..., 2]*qb[..., 1] + qa[..., 3]*qb[...,0]

    return result


def _body_frame_quaternions(times, obliquity, obl_orientation, p_rotation):
    times = np.atleast_1d(times)
    cos_obl = np.cos(obliquity)
    sin_obl = np.sqrt(1.0-cos_obl*cos_obl)
    obl = np.arccos(cos_obl)

    obl_orientation = obl_orientation
    cos_obl_orientation = np.cos(obl_orientation)
    sin_obl_orientation = np.sin(obl_orientation)

    omega_rot = 2.0*np.pi/np.exp(p_rotation)

    S = np.array([cos_obl_orientation * sin_obl, sin_obl_orientation * sin_obl, cos_obl])


    spin_rot_quat = quaternion_multiply(rotation_quaternions(np.array([0.0, 1.0, 0.0]), -obl), \
                                    rotation_quaternions(np.array([0.0, 0.0, 1.0]), -obl_orientation))

    rot_quats = rotation_quaternions(S, -omega_rot * times)

    return quaternion_multiply(spin_rot_quat, rot_quats)

def rotate_vector(rqs, v):
    nrs = rqs.shape[0]

    rqs = rqs

    vq = np.zeros((nrs, 4))
    vq[:,1:] = v

    result = quaternion_multiply(rqs, vq)
    rqs[:,1:] *= -1
    result = quaternion_multiply(result, rqs)

    return result[:,1:]

def get_visibility_illum(p_rotation = p_rotation, p_orbit=p_orbit, phi_orb = phi_orb,
                   inclination=inclination, solstice_phase=solstice_phase, obliquity=obliquity,
                   phi_rot =phi_rot, times=times):
    # A function that takes as input the factors, then returns two values: cos_insolation
    # (the parallel to the illumination) and cos_obs (the parallel to visibility)

    obl_orientation = solstice_phase
    phi_orb = phi_orb
    cos_phi_orb = np.cos(phi_orb)
    sin_phi_orb = np.sin(phi_orb)

    omega_orb = 2.0*np.pi/np.exp(p_orbit)

    no = _observer_normal_orbit_coords(inclination=inclination)
    ns = np.array([-cos_phi_orb, -sin_phi_orb, 0.0])

    orb_quats = rotation_quaternions(np.array([0.0, 0.0, 1.0]), -omega_orb*times)

    to_body_frame_quats = _body_frame_quaternions(times=times, obliquity=obliquity, obl_orientation=obl_orientation,
                                                  p_rotation=p_rotation)
    star_to_bf_quats = quaternion_multiply(to_body_frame_quats, orb_quats)

    nos = rotate_vector(to_body_frame_quats, no)
    nss = rotate_vector(star_to_bf_quats, ns)

    nside_illum = 16
    npix_illum = hp.nside2npix(nside_illum)

    pts = hp.pix2vec(nside_illum, np.arange(0, npix_illum))
    pts = np.column_stack(pts)

    cos_insolation = np.sum(nss[:, np.newaxis, :] * pts[np.newaxis, :, :], axis=2)
    cos_insolation[cos_insolation > 1] = 1.0
    cos_insolation[cos_insolation < 0] = 0.0

    cos_obs = np.sum(nos[:, np.newaxis, :] * pts[np.newaxis, :, :], axis=2)
    cos_obs[cos_obs > 1] = 1.0
    cos_obs[cos_obs < 0] = 0.0

    return cos_insolation, cos_obs

# Get the substellar and subobserver points

def get_points(p_rotation = p_rotation, p_orbit=p_orbit, phi_orb = phi_orb,
                   inclination=inclination, solstice_phase=solstice_phase, obliquity=obliquity,
                   phi_rot =phi_rot, times=times):

    obl_orientation = solstice_phase
    phi_orb = phi_orb
    cos_phi_orb = np.cos(phi_orb)
    sin_phi_orb = np.sin(phi_orb)

    omega_orb = 2.0*np.pi/np.exp(p_orbit)

    no = _observer_normal_orbit_coords(inclination=inclination)
    ns = np.array([-cos_phi_orb, -sin_phi_orb, 0.0])

    orb_quats = rotation_quaternions(np.array([0.0, 0.0, 1.0]), -omega_orb*times)

    to_body_frame_quats = _body_frame_quaternions(times=times, obliquity=obliquity, obl_orientation=obl_orientation,
                                                  p_rotation=p_rotation)
    star_to_bf_quats = quaternion_multiply(to_body_frame_quats, orb_quats)

    nos = rotate_vector(to_body_frame_quats, no)
    nss = rotate_vector(star_to_bf_quats, ns)

    return nos, nss