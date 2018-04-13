#!/bin/env python

import numpy as np
import matplotlib.pyplot as plt

import astropy.units as u
from astropy.coordinates import ICRS, Galactic

import read_mock as rm


def plot(basedir, basename, fsample):
    """
    Make a plot of the cumulative number of stars in a mock
    as a function of G band apparent magnitude.

    basedir  - directory containing the files
    basename - files are assumed to have names of the
               form <basedir>/<basename>.*.hdf5
    fsample  - use a random sampled fraction fsample
               of the stars
    """
    
    # Read the data
    data = rm.read_mock(basedir, basename, 
                        datasets=("HCoordinates",
                                  "Gmagnitude"),
                        fsample=fsample)
    hcoordinates = data["HCoordinates"]
    mag_g        = data["Gmagnitude"]
 
    # Magnitude limit depends on galactic latitude, so find latitude
    ra  = u.Quantity(hcoordinates[:,0], unit=u.radian)
    dec = u.Quantity(hcoordinates[:,1], unit=u.radian)
    coords_icrs     = ICRS(ra=ra, dec=dec)
    coords_galactic = coords_icrs.transform_to(Galactic)
    b_gal = coords_galactic.b

    plt.figure(figsize=(8.27, 11.69))

    # Make a plot of all stars with |b| < 20.0 degrees
    ind = np.abs(b_gal) <= u.Quantity(20, unit=u.degree)
    plt.subplot(2,1,1)
    plt.hist(mag_g[ind], bins=100, weights=np.ones(np.sum(ind))/fsample, cumulative=True)
    plt.yscale("log")
    plt.xlabel("Apparent magnitude in G band")
    plt.ylabel("Cumulative number of stars")
    plt.title("Stars with |b| < 20 degrees")

    # Make a plot of all stars with |b| >= 20.0 degrees
    ind = np.abs(b_gal) > u.Quantity(20, unit=u.degree)
    plt.subplot(2,1,2)
    plt.hist(mag_g[ind], bins=100, weights=np.ones(np.sum(ind))/fsample, cumulative=True)
    plt.yscale("log")
    plt.xlabel("Apparent magnitude in G band")
    plt.ylabel("Cumulative number of stars")
    plt.title("Stars with |b| > 20 degrees")
    
    plt.suptitle("Cumulative number of stars as a function of G magnitude")
    
    plt.savefig("cumulative_number.pdf")
