# Create timelines using TOAST

This is a small Python program that allows to create timelines using TOAST and
save them in a set of FITS files. It uses MPI to efficiently distribute the
workload among a set of processors.

The program requires a configuration of the instrument and of its scanning
strategy as an input, and a sky map in Healpix format. It produces a set of
FITS files, each containing a binary table HDU with as many rows as the number
of samples. The table has the following columns:

1. `TIME`, the time;
2. `THETA`, the colatitude;
3. `PHI`, the longitude;
4. `PSI`, the polarization angle;
5. `FGTOD`, the component of the signal due to foregrounds (the sky map
   provided as input);
6. `DIPTOD`, the component of the signal due to the CMB dipole;
7. `NOISETOD`, the component of the signal due to instrument noise;
8. `TOTALTOD`, the overall signal (this is *uncalibrated*);
9. `FLAGS`, the flags specific for this detector;
10. `GLFLAGS`, the global flags.

It is possible to create files for many detectors at the same time.


## Installation


### Requirements

To run this script, you must have the following Python modules available:

- TOAST, https://github.com/hpc4cmb/toast (of course!)
- NumPy, http://www.numpy.org/
- AstroPy, http://www.astropy.org/ (for FITS I/O)
- Healpy, https://github.com/healpy/healpy (for Healpix map I/O)
- A fortran compiler that is supported by `f2py` (`gfortran` is ok).


### Compilation

Download the source in a directory and run

    python setup.py build_ext --inplace

to build the Fortran extension.


## Usage

Run the program from the command line using the following syntax:

    python create_timelines.py PARAMETER_FILE OUTPUT_PATH

where `PARAMETER_FILE` is a INI file containing the parameters of the
simulation, and `OUTPUT_PATH` is the path to a directory where the FITS file
will be saved.

If you want to distribute the workload using MPI, you can use your MPI driver
(`mpirun` or whatever) as usual. In this case, the number of files saved by the
program is equal to the number of MPI processes used.


## Syntax of the parameter file

The parameter file follows the traditional syntax of INI files (see
https://docs.python.org/3/library/configparser.html).

Here is an example of a configuration file:

    [det1]
    NET = 5.23e-05
    alpha = 1.0
    fknee_Hz = 0.02
    quat = 0.0, 0.0, 1.0, 0.0
    gain_V_over_K = 35.0
    polarisation_angle_deg = 0

    [det2]
    NET = 5.23e-05
    alpha = 1.0
    fknee_Hz = 0.02
    quat = 0.0, 0.0, 1.0, 0.0
    gain_V_over_K = 35.0
    polarisation_angle_deg = 90

    [detectors]
    name_list = det1, det2

    [scanning_strategy]
    sample_rate_Hz = 84.97
    spin_period_min = 2.0
    spin_angle_deg = 65.0
    prec_period_min = 5760.0
    prec_angle_deg = 30.0
    hwp_rpm = 0.0
    hours_per_observation = 24.0
    gap_width_hours = 0.0
    chunks_per_observation = 4
    num_of_observations = 1

    [sky_model]
    sky_map = ./foreground_maps/70GHz_ecliptic.fits

The meaning of each parameter should be obvious. Keep in mind the following caveats:

1. It is possible to specify more than one detector; each of them must have its
   own section (in the example above, there are two detectors: `det1` and `det2`);
2. Only the detectors listed under the `detectors` section (see `name_list`
   above) will be actually used by the program;
3. The parameter `gain_V_over_K` is used to calibrate `TOTALTOD` from the sum
   of the columns `FGTOD`, `DIPTOD`, and `NOISETOD`. It is a float (one number)
   that is applied uniformly to all the samples in the TOD.
4. Where applicable, each parameter in tile file has the expected measure unit
   appended to the name, e.g., `prec_period_min` (the precession period) must be
   specified in minutes.
5. The parameter file used as input is saved in each FITS produced by the
   program, as a n×1 matrix of bytes.

To read the parameter file out of the primary HDU of a FITS file, you can use
something like this:

`````python
from astropy.io import fits

with fits.open('file.fits') as f:
    print(''.join([chr(x) for x in f[0].data]))
`````

## License

The code is released under a MIT license, see the file `LICENSE` for more information.
