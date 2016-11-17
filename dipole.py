#!/usr/bin/env python3
# -*- encoding: utf-8 -*-

import numpy as np
import healpy

# Colatitude of the solar system motion relative to the CMB
# (ecliptical coordinates)
SOLSYSDIR_ECL_THETA = 1.7656131194951572

# Longitude of the solar system motion relative to the CMB
# (ecliptical coordinates).
SOLSYSDIR_ECL_PHI = 2.995889600573578

# Value of the speed of the solar system motion relative to the CMB in
# m/s.
SOLSYSSPEED_M_S = 370082.2332

# Velocity vector of the Solar System motion relative to the CMB in
# m/s.
SOLSYS_SPEED_VEC_M_S = \
    SOLSYSSPEED_M_S * \
    np.array([np.sin(SOLSYSDIR_ECL_THETA) * np.cos(SOLSYSDIR_ECL_PHI),
              np.sin(SOLSYSDIR_ECL_THETA) * np.sin(SOLSYSDIR_ECL_PHI),
              np.cos(SOLSYSDIR_ECL_THETA)])


# Light speed in m/s (CODATA 2006)
SPEED_OF_LIGHT_M_S = 2.99792458e8

# Average CMB temperature in K (Mather et al. 1999, ApJ 512, 511)
T_CMB = 2.72548

################################################################################

def get_dipole_temperature(ecl_dir) -> float:
    '''Given one or more one-length versors, return the intensity of the CMB dipole

    The vectors must be expressed in the Ecliptic coordinate system.
    No kinetic dipole nor relativistic corrections are computed.
    '''
    beta = SOLSYS_SPEED_VEC_M_S / SPEED_OF_LIGHT_M_S
    gamma = (1 - np.dot(beta, beta))**(-0.5)

    return T_CMB * (1.0 / (gamma * (1 - np.dot(beta, ecl_dir))) - 1.0)

################################################################################

if __name__ == '__main__':
    nside = 32
    directions = healpy.pix2vec(nside, np.arange(healpy.nside2npix(nside)))
    healpy.write_map('dipole_map.fits', get_dipole_temperature(directions))
