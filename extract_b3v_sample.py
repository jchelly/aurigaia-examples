#!/bin/env python
#
# Find B3V stars used in figures 6 and 8 from Grand et al (2018).
# Writes out a new file in the same format as the input.
#

import numpy as np
import h5py
import read_mock as rm

import astropy.units as u
from astropy.coordinates import ICRS, Galactic


def extract_b3v_sample(basename, basedir, outfile):

    # Which quantities to read
    datasets = (
        "Vmagnitude",
        "Imagnitude",
        "VabsMagnitude", 
        "IabsMagnitude",
        "HCoordinates",
        "HCoordinatesObs", 
        "HCoordinateErrors",
        "HVelocities",
        "HVelocitiesObs",
        "HVelocityErrors",
        "Age",
        "ParticleID"
        )
    
    #
    # Selection criteria for B3V stars
    #
    # - Apparent V mag < 16.0 (filters out halo-only stars in ICC mocks)
    # - Absolute V mag in range -1.2 to -0.9
    # - V-I colour in range -0.22 to -0.16
    #
    filters = (
        rm.RangeFilter("Vmagnitude", vmax=16.0),
        rm.RangeFilter("VabsMagnitude", vmin=-1.2, vmax=-0.9),
        rm.DifferenceFilter("VabsMagnitude","IabsMagnitude", vmin=-0.22, vmax=-0.16)
        )

    # Read the mock
    data = rm.read_mock(basedir, basename, datasets=datasets, filters=filters)
    nstars = data["Vmagnitude"].shape[0]

    # Write the output file
    with h5py.File(outfile, "w") as outfile:

        # Write data arrays
        grp = outfile.create_group("Stardata")
        for name in datasets:
            grp[name] = data[name]

        # Update header with new number of stars etc
        data["Header"]["NumFilesPerSnapshot"] = 1
        if "NumPart_ThisFile" in data["Header"]:
            data["Header"]["NumPart_ThisFile"] = nstars
            data["Header"]["NumPart_Total"]    = nstars
        if "NumStars_ThisFile" in data["Header"]:
            data["Header"]["NumStars_ThisFile"] = nstars
            data["Header"]["NumStars_Total"]    = nstars        

        # Write header and parameters
        grp = outfile.create_group("Header")
        for name in data["Header"]:
            grp.attrs[name] = data["Header"][name]
        grp = outfile.create_group("Parameters")
        for name in data["Parameters"]:
            grp.attrs[name] = data["Parameters"][name]




if __name__ == "__main__":

    basedir  = "/cosma5/data/jch/Gaia/HITS/v1/kroupaIMF/level3_MHD/halo_6/mockdir_angle030/"
    basename = "mock_030"
    outfile  = "H_Au06_B3V_Ex.hdf5"

    #basedir  = "/cosma5/data/jch/Gaia/HITS/v1/kroupaIMF/level3_MHD/halo_6/mockdir_angle030/"
    #basename = "mock_noex_030"
    #outfile  = "H_Au06_B3V_NoEx.hdf5"

    #basedir  = "/cosma5/data/jch/Gaia/ICC/v1/chabrierIMF/level3_MHD/halo_6/mockdir_angle030/"
    #basename = "mock_noex_030"
    #outfile  = "I_Au06_B3V.hdf5"

    extract_b3v_sample(basename, basedir, outfile)
