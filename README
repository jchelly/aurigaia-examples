# AuriGaia Examples for Python

These python scripts demonstrate how to read the AuriGaia mock catalogue
files and make some basic plots. We also include a brief description of
the file format here.


## Dependencies

The following modules are required to run these scripts:

  * numpy - for efficient handling of large arrays
  * h5py - for reading the HDF5 files
  * astropy - for transformations between coordinate systems
  * healpy - for making healpix maps
  * matplotlib - for making plots

## Mock file format

### Layout of the files

Each mock catalogue consists of a set of 96 files. They have names of the
form

  mock_XXX.Y.hdf5 (HITS mocks, with dust extinction) or
  mock_noex_XXX.Y.hdf5 (ICC mocks, without dust extinction)

where XXX is the angle between the major axis of the bar of the galaxy and
the Sun-Galactic centre line, and Y is the file index which runs from 0 to 95.

The files are in HDF5 format. Each file contains a Stardata group which
contains a number of datasets with the properties of the synthetic stars.
Each dataset has attributes "Description" and "Units" which contain a 
description of the meaning of the dataset and its units, respectively.


### True vs observed quantities

  * Datasets with the suffix "Obs" have been randomly perturbed in line with
    the magnitude of the errors expected in Gaia DR2.

  * For each dataset with the "Obs" suffix there's a corresponding dataset 
    with the suffix "Error", which contains the magnitude of the expected 
    error.

  * For each "Obs" dataset there is also a dataset without the "Obs" suffix.
    This contains the "true" value before errors are applied.


### Positions of the stars

The datasets HCoordinatesObs and HCoordinates contain the equatorial
coordinates of the stars with and without errors, respectively. The
contents of these datasets are as follows:

  HCoordinates[:,0] - right ascension in radians
  HCoordinates[:,1] - declination in radians
  HCoordinates[:,2] - parallax in arcseconds

HCoordinatesObs contains the corresponding values with errors, and 
HCoordinateErrors contains the errors on these quantities.

cartesian_coords.py contains a function which shows how to convert these 
positions into cartesian coordinates using astropy.


### Velocities of the stars

The datasets HVelocitiesObs and HVelocities contain the proper motions
and radial velocities of the stars with and without errors. The contents
of these datasets are as follows:

  HVelocities[:,0] - proper motion in right ascension*cos(declination), 
                     in radians
  HVelocities[:,1] - declination in radians
  HVelocities[:,2] - radial velocity in km/sec

HVelocitiesObs contains the corresponding values with errors, and 
HVelocityErrors contains the errors on these quantities.

cartesian_coords.py contains a function which shows how to convert these 
velocities into Galactic cartesian coordinates using astropy.


### Magnitudes

Apparent magnitudes in the Gaia G, GRP and GBP bands are stored in their
own datasets with and without errors according to the naming scheme
described above.

The mocks also contain apparent magnitudes in the UBVRIJHK bands in the
dataset Magnitudes. This is a 2D array defined as follows:

Magnitudes[:,0] - apparent magnitude in U
Magnitudes[:,1] - apparent magnitude in B
Magnitudes[:,2] - apparent magnitude in R
Magnitudes[:,3] - apparent magnitude in J
Magnitudes[:,4] - apparent magnitude in H
Magnitudes[:,5] - apparent magnitude in K
Magnitudes[:,6] - apparent magnitude in V
Magnitudes[:,7] - apparent magnitude in I


## Reading the catalogues

### Using h5py directly

Individual files can be read as follows:

```
  import h5py

  f = h5py.File("mock_noex_030.58.hdf5","r")
  age = f["Stardata/Age"][...]
  mag_g = f["Stardata/Gmagnitude"][...]
  hcoordinates = f["Stardata/HCoordinates"][...]
  ...
  f.close()
```

Using the [...] notation ensures that h5py reads the whole dataset into
memory and returns a numpy array rather than a h5py.Dataset object.


### Using the read_aurigaia() function

The partitioning of stars between files is arbitrary, so in general it
will be necessary to read all of the files in a mock catalogue to make use 
of the data. The file read_aurigaia.py provides a function read_aurigaia() to do 
this while only reading quantities which have been specifically asked for 
and optionally filtering out unwanted stars as it goes to save memory.

This can be used as follows:

```
import read_aurigaia as ra

data = ra.read_aurigaia(basedir="./ICC/v2/chabrierIMF/level3_MHD/halo_6/mockdir_angle030",
                        basename="mock_noex_030",
                        datasets=("Age","Gmagnitude","HCoordinates","HVelocities"))
```

This will return a dictionary 'data' where the keys are the dataset names
specified in the datasets parameter and the values are numpy arrays with
the data read from all of the mock files.

It's also possible to filter out unwanted stars using the 'filters' parameter,
which specifies a list of classes to use to filter the stars from each file
as it is read.

For example, to extract a 1% sample of stars from a set of mock files:

```
data = ra.read_aurigaia(basedir="./ICC/v2/chabrierIMF/level3_MHD/halo_6/mockdir_angle030",
                        basename="mock_noex_030",
                        datasets=("Age","Gmagnitude","HCoordinates","HVelocities"),
                        filters=(ra.RandomSampleFilter(0.01)))
```

See read_aurigaia.py for examples of how to define these filters. The example
extract_b3v_sample.py shows how to read in the B3V stars from figures 6 and
8 of Grand et al by applying several filters in sequence.


## Example scripts

### Making a plot in cartesian coordinates

The script cartesian_plot.py reads in a random sample of the stars from a
mock, converts to (heliocentric) galactic cartesian coordinates, and makes 
a plot of the positions of the stars.

The script takes two command line parameters: the directory containing the
mock files and the base name of the files. E.g:
```
python ./cartesian_plot.py ./ICC/v2/chabrierIMF/level3_MHD/halo_6/mockdir_angle030 mock_noex_030
```
The resulting plot is saved to a .png file.


### Making a colour magnitude diagram

The script colour_magnitude.py reads in a random sample of the stars from a
mock and plots a colour magnitude diagram.

The script takes two command line parameters: the directory containing the
mock files and the base name of the files. E.g:
```
python ./colour_magnitude.py ./ICC/v2/chabrierIMF/level3_MHD/halo_6/mockdir_angle030 mock_noex_030
```
The resulting plot is saved to a .pdf file.


### Making a HEALPix map

The example healpix_age.py shows how to make a HEALPix map from the mocks.
Each pixel is coloured by the mean age of the stars in that pixel. The
functions healpy.pixelfunc.ang2pix and numpy.bincount can be used to
efficiently map stars to HEALPix pixels.

As before, this can be run with
```
python ./healpix_age.py ./ICC/v2/chabrierIMF/level3_MHD/halo_6/mockdir_angle030 mock_noex_030
```
and the result is written to a .pdf file.


### Extracting a sample of stars

The example extract_b3v_sample.py shows how to use the filters parameter
of read_aurigaia() to select the sample of B3V stars shown in figures 6
and 8 of Grand et al.

The script can be run with
```
python ./extract_b3v_sample.py ./ICC/v2/chabrierIMF/level3_MHD/halo_6/mockdir_angle030 mock_noex_030 sample_B3V.hdf5
```
where the third command line parameter is the name of the output file to
create. This file is in the same format as the input mock files.
