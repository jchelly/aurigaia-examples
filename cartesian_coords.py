#!/bin/env python

import numpy as np
import h5py

import astropy.units as u
from astropy.coordinates import ICRS, Galactic


def galactic_cartesian_coordinates(hcoordinates, hvelocities=None):
    """
    Calculate Galactic cartesian coordinates of stars
    from a mock catalogue given equatorial (ICRS) coordinates.
    In this system the disk of the galaxy is in the x-y
    plane with the Sun at the origin and the galactic
    centre in the +x direction.

    hcoordinates - [N,3] array of ICRS equatorial coordinates of N 
                   stars  Coordinates are defined as follows:

    hcoordinates[:,0] - right acension of the stars in radians
    hcoordinates[:,1] - declination of the stars in radians
    hcoordinates[:,2] - parallax of the stars in arcseconds

    Can optionally also specify velocities:

    hvelocities[:,0] - proper motion in right ascension * cos(declination) 
                       (arcsec/yr)
    hvelocities[:,1] - proper motion in declination (arcsec/yr)
    hvelocities[:,2] - radial velocity (km/s)

    Returns:

    pos - [N,3] array with heliocentric cartesian coordinates in kpc.

    If hvelocities is not None, also returns

    vel - [N,3] array with heliocentric cartesian velocities in km/s
    """

    # Extract coordinate components from the HCoordinates dataset
    ra       = u.Quantity(hcoordinates[:,0], unit=u.radian)
    dec      = u.Quantity(hcoordinates[:,1], unit=u.radian)
    parallax = u.Quantity(hcoordinates[:,2], unit=u.arcsec)

    # Extract velocities from HVelocities dataset if present
    if hvelocities is not None:
        pm_ra_cosdec = u.Quantity(hvelocities[:,0], unit=u.arcsec/u.year)
        pm_dec       = u.Quantity(hvelocities[:,1], unit=u.arcsec/u.year)
        rv           = u.Quantity(hvelocities[:,2], unit=u.km/u.second)

    # Calculate distance to each star:
    # Distance in parsecs is just 1.0/(parallax in arcsecs),
    # but here we let astropy deal with the units.
    dist = parallax.to(u.kpc, equivalencies=u.parallax())

    # Translate to galactic coordinates using astropy
    if hvelocities is not None:
        # Have positions and velocities
        coords = ICRS(ra=ra, dec=dec, distance=dist,
                      pm_ra_cosdec=pm_ra_cosdec, pm_dec=pm_dec, 
                      radial_velocity=rv).transform_to(Galactic)
        pos = coords.cartesian.get_xyz().transpose()
        vel = coords.cartesian.differentials["s"].get_d_xyz().transpose()
        return pos, vel
    else:
        # Just have positions
        coords = ICRS(ra=ra, dec=dec, distance=dist).transform_to(Galactic)
        pos = coords.cartesian.get_xyz().transpose()
        return pos
