#!/bin/env python

from __future__ import print_function

import numpy as np
import h5py


def read_mock(basedir, basename, datasets, fsample=1.0):
    """
    Read in the specified quantities from a set of mock
    files.

    basedir  - directory containing the files
    basename - files are assumed to have names of the
               form <basedir>/<basename>.*.hdf5
    datasets - a sequence of strings with the names of
               the HDF5 datasets to read from the
               Stardata groups in the mock files.
    fsample  - return a random sampled fraction fsample
               of the stars

    Returns a dictionary where the keys are dataset names
    and the values are numpy arrays containing the data.
    This dictionary also contains "Header" and "Parameters"
    entries which contain the values from the Header and
    Parameters groups in the HDF5 files.
    """

    # Create the output dictionary
    data = {name : [] for name in datasets}
    
    # Loop over input files
    ifile  = 0
    nfiles = 1
    while ifile < nfiles:

        # Open the next file
        fname = "%s/%s.%d.hdf5" % (basedir, basename, ifile)
        print("Reading: ", fname)
        with h5py.File(fname, "r") as infile:
            # Get number of files from first file
            if ifile == 0:
                nfiles = infile["Header"].attrs["NumFilesPerSnapshot"]
                # Copy header to output
                data["Header"] = {}
                for name in infile["Header"].attrs:
                    data["Header"][name] = infile["Header"].attrs[name]
                # Copy parameters to output
                data["Parameters"] = {}
                for name in infile["Parameters"].attrs:
                    data["Parameters"][name] = infile["Parameters"].attrs[name]
            # Read datasets
            for name in datasets:
                if fsample >= 1.0:
                    data[name].append(infile["Stardata"][name][...])
                else:
                    ind = np.random.rand(infile["Stardata"][name].shape[0]) < fsample
                    data[name].append(infile["Stardata"][name][...][ind,...])

        # Next file
        ifile += 1

    # Combine arrays
    for name in datasets:
        data[name] = np.concatenate(data[name], axis=0)
    
    return data
