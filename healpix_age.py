#!/bin/env python

import matplotlib.pyplot as plt
import numpy as np
import healpy as hp

import read_mock as rm


def plot(basedir, basename, fsample=0.1):
    """
    Make a healpix map of the age of stars in a mock.

    basedir  - directory containing the files
    basename - files are assumed to have names of the
               form <basedir>/<basename>.*.hdf5
    fsample  - use a random sampled fraction fsample
               of the stars
    """
    
    # Read the data
    data = rm.read_mock(basedir, basename, 
                        datasets=("HCoordinates",
                                  "Age"),
                        filters=(rm.RandomSampleFilter(fsample),))

    ra_radians  = data["HCoordinates"][:,0]
    dec_radians = data["HCoordinates"][:,1]
    age_gyr     = data["Age"]
     
    # Find healpix pixel index for each star
    NSIDE = 512
    ipix = hp.pixelfunc.ang2pix(NSIDE, np.degrees(ra_radians), np.degrees(dec_radians), lonlat=True)

    # Find mean age of stars in each pixel
    mean_age = (np.bincount(ipix, minlength=hp.nside2npix(NSIDE), weights=age_gyr) / 
                np.bincount(ipix, minlength=hp.nside2npix(NSIDE)))

    plt.figure(figsize=(8.27, 11.69))

    # Make a plot in equatorial coordinates
    plt.subplot(2,1,1)
    hp.mollview(mean_age, coord="C", fig=plt.gcf().number, hold=True)
    hp.visufunc.graticule()

    # Make a plot in galactic coordinates
    plt.subplot(2,1,2)
    hp.mollview(mean_age, coord="CG", fig=plt.gcf().number, hold=True)
    hp.visufunc.graticule()

    plt.savefig("healpix_age.pdf")
    
