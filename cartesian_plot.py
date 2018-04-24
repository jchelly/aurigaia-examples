#!/bin/env python

import matplotlib.pyplot as plt
import matplotlib.colors as colors
import matplotlib.cm as cm

import numpy as np

import read_mock as rm
import cartesian_coords as cc


def plot(basedir, basename, fsample=0.01):
    """
    Make a plot of the specified mock in cartesian
    coordinates showing number density of stars.
    
    basedir  - directory containing the files
    basename - files are assumed to have names of the
               form <basedir>/<basename>.*.hdf5
    fsample  - use a random sampled fraction fsample
               of the stars

    NOTE: the transformation to cartesian coordinates is
    not reliable when there are significant errors in the
    parallax, so we use the true coordinates here.
    """

    # Read coordinates from the mock
    data = rm.read_mock(basedir, basename, 
                        datasets=("HCoordinates",),
                        fsample=fsample)
    hcoordinates = data["HCoordinates"]
 
    # Get position of the Sun in kpc relative to the galactic centre
    sun_pos = np.asarray((data["Parameters"]["solarradius"],
                          data["Parameters"]["solarheight"],
                          0.0), dtype=float) * 1000.0
    # Position of the Galactic centre relative to the Sun
    gc_pos = -sun_pos

    # Transform into cartesian coordinates
    pos = cc.galactic_cartesian_coords(hcoordinates)
 
    # Make the plot
    plt.figure(figsize=(8.27, 11.69))

    # x-y plot
    plt.subplot(2,2,1).set_aspect("equal")
    plt.hist2d(pos[:,0], pos[:,1], range=((-20,20),(-20,20)), bins=(500,500), 
               norm=colors.LogNorm(1,1000), cmap=cm.Greys)
    plt.xlabel("x [kpc]")
    plt.ylabel("y [kpc]")
    plt.plot(0,         0,         "r.", label="Sun")
    plt.plot(gc_pos[0], gc_pos[1], "g+", label="Galactic centre")

    # x-z plot
    plt.subplot(2,2,2).set_aspect("equal")
    plt.hist2d(pos[:,0], pos[:,2], range=((-20,20),(-20,20)), bins=(500,500), 
               norm=colors.LogNorm(1,1000), cmap=cm.Greys)
    plt.xlabel("x [kpc]")
    plt.ylabel("y [kpc]")
    plt.plot(0,         0,         "r.", label="Sun")
    plt.plot(gc_pos[0], gc_pos[2], "g+", label="Galactic centre")

    # y-z plot
    plt.subplot(2,2,3).set_aspect("equal")
    plt.hist2d(pos[:,1], pos[:,2], range=((-20,20),(-20,20)), bins=(500,500), 
               norm=colors.LogNorm(1,1000), cmap=cm.Greys)
    plt.xlabel("y [kpc]")
    plt.ylabel("z [kpc]")
    plt.plot(0,         0,         "r.", label="Sun")
    plt.plot(gc_pos[1], gc_pos[2], "g+", label="Galactic centre")
    plt.legend(loc="lower left")

    plt.suptitle("Projected number density of stars")

    #plt.savefig("cartesian_plot.pdf")
    plt.savefig("cartesian_plot.png", dpi=150)
