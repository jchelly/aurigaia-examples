#!/bin/env python

import numpy as np
import matplotlib.pyplot as plt
import matplotlib.colors as colors

import read_aurigaia as ra


def plot(basedir, basename, fsample):
    """
    Plot the color magnitude diagram for all stars
    """
    
    # Read the data
    data = ra.read_aurigaia(basedir, basename, 
                            datasets=("VabsMagnitude",
                                      "IabsMagnitude",
                                      "Magnitudes"),
                            filters=(ra.RandomSampleFilter(fsample),))

    mag_v_abs = data["VabsMagnitude"]
    mag_i_abs = data["IabsMagnitude"]
    mag_v_app = data["Magnitudes"][:,6]

    # Select only stars with V<16
    ind = mag_v_app < 16.0

    # Make a 2D histogram
    plt.hist2d(mag_v_abs[ind]-mag_i_abs[ind], mag_v_abs[ind], bins=200, range=((-1,4),(-5,10)), norm=colors.LogNorm())
    plt.ylim(10,-5)
    plt.xlabel("V-I")
    plt.ylabel("V")

    plt.savefig("cmd.pdf")


if __name__ == "__main__":

    import sys

    basedir  = sys.argv[1]
    basename = sys.argv[2]

    plot(basedir, basename, fsample=0.1)
