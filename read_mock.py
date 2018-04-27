#!/bin/env python

import numpy as np
import h5py


class BaseFilter():
    """
    Base class for filters for use with read_mock.
    This just provides the loop over datasets.
    Implementations need to provide __init__ and
    stars_to_keep.
    """
    def __call__(self, data):
        keep = None
        for name in data:
            if name not in ("Header", "Parameters"):
                if keep is None:
                    nstars = data[name].shape[0]
                    keep = self.stars_to_keep(nstars, data)
                data[name] = data[name][keep,...]
        return data
    def _in_range(self, arr, vmin, vmax):
        if vmin is not None and vmax is not None:
            return np.logical_and(arr>=self.vmin, arr<self.vmax)
        elif vmin is not None:
            return arr>=self.vmin
        elif vmax is not None:
            return arr<self.vmax
        else:
            raise ValueError("Must specify vmin or vmax!")


class RandomSampleFilter(BaseFilter):
    """
    This is intended for use with read_mock() - 
    set filter=RandomSampleFilter(fsample) to get
    a random sample of the stars.
    """
    def __init__(self, fsample):
        self.fsample = fsample
    def stars_to_keep(self, nstars, data):
        return np.random.rand(nstars) < self.fsample


class RangeFilter(BaseFilter):
    """
    This is intended for use with read_mock() - 
    set filter=RangeFilter(quantity,vmin,vmax) to get
    only stars with quantity 'quantity' in range
    vmin to vmax.
    """
    def __init__(self, quantity, vmin=None, vmax=None):
        self.quantity = quantity
        self.vmin = vmin
        self.vmax = vmax
    def stars_to_keep(self, nstars, data):
        return self._in_range(data[self.quantity], self.vmin, self.vmax)


class DifferenceFilter(BaseFilter):
    """
    This is intended for use with read_mock() - 
    set filter=DifferenceFilter(quantity1,quantity2,vmin,vmax)
    to get only stars with quantitiy1-quantity2 in range
    vmin to vmax.
    """
    def __init__(self, quantity1, quantity2, vmin=None, vmax=None):
        self.quantity1 = quantity1
        self.quantity2 = quantity2
        self.vmin = vmin
        self.vmax = vmax
    def stars_to_keep(self, nstars, data):
        diff = data[self.quantity1]-data[self.quantity2]
        return self._in_range(diff, self.vmin, self.vmax)


def read_mock(basedir, basename, datasets, filters=None):
    """
    Read in the specified quantities from a set of mock
    files.

    basedir  - directory containing the files
    basename - files are assumed to have names of the
               form <basedir>/<basename>.*.hdf5
    datasets - a sequence of strings with the names of
               the HDF5 datasets to read from the
               Stardata groups in the mock files.
    filters  - each callable in this sequence is called on 
               the data from each file. It should return a copy
               of the data dict with unwanted stars filtered out.

    Returns a dictionary where the keys are dataset names
    and the values are numpy arrays containing the data.
    This dictionary also contains "Header" and "Parameters"
    entries which contain the values from the Header and
    Parameters groups in the HDF5 files.
    """

    # Synthetic dataset names - these will be extracted from
    # the Magnitudes dataset if necessary
    extra_mag_index = {}
    for iband, band in enumerate(("U","B","R","J","H","K","V","I")):
        extra_mag_index[band+"magnitude"] = iband

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
            # Read data from this file
            mag_data = None
            data_thisfile = {}
            for name in datasets:
                if name in extra_mag_index:
                    # This is a magnitude stored in the Magnitudes dataset
                    if mag_data is None:
                        mag_data = infile["Stardata"]["Magnitudes"][...]
                    data_thisfile[name] = mag_data[:,extra_mag_index[name]]
                else:
                    # This is a normal dataset
                    data_thisfile[name] = infile["Stardata"][name][...]
            data_thisfile["Header"]     = data["Header"]
            data_thisfile["Parameters"] = data["Parameters"]
            # Apply filter function(s), if any
            if filters is not None:
                for filter in filters:
                    data_thisfile = filter(data_thisfile)
            # Append data to output
            for name in datasets:
                data[name].append(data_thisfile[name])

        # Next file
        ifile += 1

    # Combine arrays
    for name in datasets:
        data[name] = np.concatenate(data[name], axis=0)
    
    return data
