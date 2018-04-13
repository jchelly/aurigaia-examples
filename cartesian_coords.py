#!/bin/env python

import numpy as np
import h5py

import astropy.units as u
from astropy.coordinates import ICRS, Galactic


def galactic_cartesian_coords(hcoordinates):
    """
    Calculate Galactic cartesian coordinates of stars
    from a mock catalogue given equatorial coordinates.
    In this system the disk of the galaxy is in the x-y
    plane with the Sun at the origin and the galactic
    centre in the +x direction.

    hcoordinates - [N,3] array of ICRS equatorial coordinates of N 
                   stars  Coordinates are defined as follows:

    hcoordinates[:,0] - right acension of the stars in radians
    hcoordinates[:,1] - declination of the stars in radians
    hcoordinates[:,2] - parallax of the stars in arcseconds

    Returns:

    pos - [N,3] array with heliocentric cartesian coordinates in kpc.
    """

    # Extract coordinate components from the HCoordinates dataset
    ra       = u.Quantity(hcoordinates[:,0], unit=u.radian)
    dec      = u.Quantity(hcoordinates[:,1], unit=u.radian)
    parallax = u.Quantity(hcoordinates[:,2], unit=u.arcsec)

    # Calculate distance to each star:
    # Distance in parsecs is just 1.0/(parallax in arcsecs),
    # but here we let astropy deal with the units.
    dist = parallax.to(u.kpc, equivalencies=u.parallax())

    # Translate ra and dec into galactic coordinates using astropy
    # (this is just a rotation)
    coords_icrs     = ICRS(ra=ra, dec=dec, distance=dist)
    coords_galactic = coords_icrs.transform_to(Galactic)
    
    # Convert galactic lattitude/longitude/distance to cartesian coordinates
    return coords_galactic.cartesian.get_xyz().transpose()
