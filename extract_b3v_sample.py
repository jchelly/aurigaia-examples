#!/bin/env python
#
# Find B3V stars used in figures 6 and 8 from Grand et al (2018).
# Writes out a new file in the same format as the input.
#

import numpy as np
import h5py
import read_aurigaia as ra


def extract_b3v_sample(basedir, basename, outfile):

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
        ra.RangeFilter("Vmagnitude", vmax=16.0),
        ra.RangeFilter("VabsMagnitude", vmin=-1.2, vmax=-0.9),
        ra.DifferenceFilter("VabsMagnitude","IabsMagnitude", vmin=-0.22, vmax=-0.16)
        )

    # Read the mock
    data = ra.read_aurigaia(basedir, basename, datasets=datasets, filters=filters)
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

    import sys
    
    basedir  = sys.argv[1]
    basename = sys.argv[2]
    outfile  = sys.argv[3]

    extract_b3v_sample(basedir, basename, outfile)
