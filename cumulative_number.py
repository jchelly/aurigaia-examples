#!/bin/env python

import numpy as np
import matplotlib.pyplot as plt

import astropy.units as u
from astropy.coordinates import ICRS, Galactic

import read_mock as rm


def plot(basedir, basename, fsample=0.1):
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
                        filters=(rm.RandomSampleFilter(fsample),))

    hcoordinates = data["HCoordinates"]
    mag_g        = data["Gmagnitude"]
 
    # Get equatorial coordinates
    ra  = u.Quantity(hcoordinates[:,0], unit=u.radian)
    dec = u.Quantity(hcoordinates[:,1], unit=u.radian)
    parallax = u.Quantity(hcoordinates[:,2], unit=u.arcsec)

    # Calculate distance to each star:
    # Distance in parsecs is just 1.0/(parallax in arcsecs),
    # but here we let astropy deal with the units.
    dist = parallax.to(u.kpc, equivalencies=u.parallax())
    
    # Calculate galactic latitude of each star
    coords = ICRS(ra=ra, dec=dec, distance=dist).transform_to(Galactic)
    b_gal = coords.b

    plt.figure(figsize=(8.27, 11.69))

    # Make a plot of all stars with |b| < 20.0 degrees
    ind = np.abs(b_gal) <= u.Quantity(20, unit=u.degree)
    plt.subplot(2,2,1)
    plt.hist(mag_g[ind], bins=100, range=(0,21), weights=np.ones(np.sum(ind))/fsample, cumulative=True, histtype="step")
    plt.yscale("log")
    plt.xlabel("Apparent magnitude in G band")
    plt.ylabel("Cumulative number of stars")
    plt.xlim(0,21)
    plt.ylim(1.0e0,1.0e9)
    plt.title("|b| < 20 degrees")

    # Make a plot of all stars with |b| >= 20.0 degrees
    ind = np.abs(b_gal) > u.Quantity(20, unit=u.degree)
    plt.subplot(2,2,2)
    plt.hist(mag_g[ind], bins=100, range=(0,21), weights=np.ones(np.sum(ind))/fsample, cumulative=True, histtype="step")
    plt.yscale("log")
    plt.xlabel("Apparent magnitude in G band")
    plt.ylabel("Cumulative number of stars")
    plt.xlim(0,21)
    plt.ylim(1.0e0,1.0e9)
    plt.title("|b| > 20 degrees")
    
    # Make a plot of all stars with |b| < 20.0 degrees
    ind = np.abs(b_gal) <= u.Quantity(20, unit=u.degree)
    plt.subplot(2,2,3)
    plt.hist(dist[ind], bins=10**np.linspace(-3.0,3.0,100), weights=np.ones(np.sum(ind))/fsample, cumulative=True, histtype="step")
    plt.xscale("log")
    plt.yscale("log")
    plt.xlabel("Distance [kpc]")
    plt.ylabel("Cumulative number of stars")
    plt.xlim(1.0e-3,1.0e3)
    plt.ylim(1.0,1.0e9)
    plt.title("|b| < 20 degrees")

    # Make a plot of all stars with |b| > 20.0 degrees
    ind = np.abs(b_gal) > u.Quantity(20, unit=u.degree)
    plt.subplot(2,2,4)
    plt.hist(dist[ind], bins=10**np.linspace(-3.0,3.0,100), weights=np.ones(np.sum(ind))/fsample, cumulative=True, histtype="step")
    plt.xscale("log")
    plt.yscale("log")
    plt.xlabel("Distance [kpc]")
    plt.ylabel("Cumulative number of stars")
    plt.xlim(1.0e-3,1.0e3)
    plt.ylim(1.0,1.0e9)
    plt.title("|b| > 20 degrees")


    plt.suptitle("Cumulative number of stars as a function of G magnitude")
    
    plt.savefig("cumulative_number.pdf")
