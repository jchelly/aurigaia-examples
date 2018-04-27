#!/bin/env python

import matplotlib.pyplot as plt
import numpy as np
import healpy as hp

import read_mock as rm


def plot(basedir, basename, fsample=0.1):
    """
    Make a healpix map of the number density of stars in
    a mock.

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

    ra_radians  = data["HCoordinates"][:,0]
    dec_radians = data["HCoordinates"][:,1]
    mag_g       = data["Gmagnitude"]
     
    # Find healpix pixel index for each star
    NSIDE = 512
    ipix = hp.pixelfunc.ang2pix(NSIDE, np.degrees(ra_radians), np.degrees(dec_radians), lonlat=True)

    # Count stars in each healpix pixel
    num_per_pix = np.bincount(ipix, minlength=hp.nside2npix(NSIDE)).astype(float)/fsample

    plt.figure(figsize=(8.27, 11.69))

    # Make a plot in equatorial coordinates
    plt.subplot(2,1,1)
    hp.mollview(np.log10(num_per_pix+1), coord="C", fig=plt.gcf().number, hold=True, min=0.5, max=3.7)
    hp.visufunc.graticule()

    # Make a plot in galactic coordinates
    plt.subplot(2,1,2)
    hp.mollview(np.log10(num_per_pix+1), coord="CG", fig=plt.gcf().number, hold=True, min=0.5, max=3.7)
    hp.visufunc.graticule()

    plt.savefig("healpix_number_density.pdf")
    
