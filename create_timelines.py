#!/usr/bin/env python3
# -*- encoding: utf-8 -*-

'''Usage: {basename} PARAMETER_FILE OUTPUT_PATH
'''

from collections import namedtuple
from configparser import ConfigParser
from typing import Any, Dict, List, Union
import io
import logging as log
import os
import sys

import astropy.io.fits as fits
import healpy
import numpy as np  # typing: ignore
import dipole as dip
import toast
import toast.tod as tt
import toast.tod.qarray as qa

__version__ = '0.3.0'

# Constantly slewing precession axis
DEG_ANGLE_PER_DAY = 360.0 / 365.25

##############################################################################

# Information about one CORE detector
Detector = namedtuple('Detector',
                      ['name',
                       'net',
                       'alpha',
                       'fknee_Hz',
                       'polarisation_angle_deg',
                       'gain',
                       'quat'])

# Information about the simulation as a whole
Configuration = namedtuple('Configuration',
                           ['detectors',
                            'sample_rate_Hz',
                            'spin_period_min',
                            'spin_angle_deg',
                            'prec_period_min',
                            'prec_angle_deg',
                            'hwp_rpm',
                            'hours_per_obs',
                            'gap_width_hours',
                            'chunks_per_obs',
                            'num_obs',
                            'nside',
                            'sky_map',
                            'sky_map_hdu',
                            'i_column',
                            'q_column',
                            'u_column',
                            'param_file_contents'])

################################################################################

def parse_vector(s: str) -> List[float]:
    'Transform a vector in the form "1,2,3" into a list of integers [1, 2, 3]'
    return np.array([float(x) for x in s.split(',')])

################################################################################

def int_or_str(s: str) -> Union[int, str]:
    '''Try to convert "s" into an integer. If the conversion fails, return "s".'''
    try:
        return int(s)
    except ValueError:
        return s

################################################################################

def build_configuration_obj(conf_file: ConfigParser) -> Configuration:
    'Return a "Configuration" object built using the values in "conf_file"'
    try:
        # Step 1: determine how many detectors are defined in the INI file
        det_names = sorted([x.strip()
                            for x in conf_file['detectors']['name_list'].split(',')])

        # Step 2: build a "Detector" structure for each of them, taking
        # the parameters from the corresponding section in the parameter
        # file
        detectors = [] # type: List[Detector]
        for name in det_names:
            cur_det = conf_file[name]
            detectors.append(Detector(name=name,
                                      net=cur_det.getfloat('NET'),
                                      polarisation_angle_deg=cur_det.getfloat('polarisation_angle_deg', 0.0),
                                      alpha=cur_det.getfloat('alpha'),
                                      fknee_Hz=cur_det.getfloat('fknee_Hz'),
                                      gain=cur_det.getfloat('gain_V_over_K'),
                                      quat=parse_vector(cur_det.get('quat'))))

        # Step 3: build the global "Configuration" object
        scanning = conf_file['scanning_strategy']
        sky_model = conf_file['sky_model']

        # Save a copy of the parameter file into the Configuration object
        param_file_contents = io.StringIO()
        conf_file.write(param_file_contents)
        param_file_contents = np.array(
            list(param_file_contents.getvalue().encode('utf-8')))

        result = Configuration(detectors=detectors,
                               sample_rate_Hz=scanning.getfloat('sample_rate_Hz'),
                               spin_period_min=scanning.getfloat('spin_period_min'),
                               spin_angle_deg=scanning.getfloat('spin_angle_deg'),
                               prec_period_min=scanning.getfloat('prec_period_min'),
                               prec_angle_deg=scanning.getfloat('prec_angle_deg'),
                               hwp_rpm=scanning.getfloat('hwp_rpm'),
                               hours_per_obs=scanning.getfloat('hours_per_observation'),
                               gap_width_hours=scanning.getfloat('gap_width_hours'),
                               chunks_per_obs=scanning.getint('chunks_per_observation'),
                               num_obs=scanning.getint('num_of_observations'),
                               nside=sky_model.getint('nside', 4096),
                               sky_map=sky_model.get('sky_map'),
                               sky_map_hdu=int_or_str(sky_model.get('map_hdu', 1)),
                               i_stokes=int_or_str(sky_model.get('i_stokes', 0)),
                               q_stokes=int_or_str(sky_model.get('q_stokes', 1)),
                               u_stokes=int_or_str(sky_model.get('u_stokes', 2)),
                               param_file_contents=param_file_contents)
    except KeyError as err:
        log.error('section/key not found: %s', err)
        sys.exit(1)

    return result

################################################################################

def main(parameter_file, outdir):
    'The main function'

    # Set up the logging system
    logformat = '[%(asctime)s %(levelname)s] %(message)s'
    log.basicConfig(level=log.INFO,
                    format=logformat)

    # This is the 2-level toast communicator.  By default,
    # there is just one group which spans MPI_COMM_WORLD.
    comm = toast.Comm()

    if comm.comm_world.rank == 0:
        log.info("running with %d processes", comm.comm_world.size)

    # Read the parameter file
    conf = None
    if comm.comm_world.rank == 0:
        configuration_parser = ConfigParser()
        configuration_parser.read(parameter_file)
        conf = build_configuration_obj(configuration_parser)
        log.info('configuration read from file "%s"', parameter_file)
    conf = comm.comm_world.bcast(conf, root=0)

    if outdir != '.':
        if not os.path.isdir(outdir):
            os.mkdir(outdir)

    # construct the list of intervals which split the whole observation time
    intervals = tt.regular_intervals(n=int(conf.num_obs),
                                     start=0.0,
                                     first=0,
                                     rate=conf.sample_rate_Hz,
                                     duration=3600*float(conf.hours_per_obs),
                                     gap=3600*float(conf.gap_width_hours))

    detquats = {}  # type: Dict[str, Any]
    for det in conf.detectors:
        detquats[det.name] = det.quat

    # Create the noise model for this observation

    rate = {}  # type: Dict[str, float]
    fmin = {}  # type: Dict[str, float]
    fknee = {}  # type: Dict[str, float]
    alpha = {}  # type: Dict[str, float]
    net = {}  # type: Dict[str, float]
    for det in conf.detectors:
        rate[det.name] = conf.sample_rate_Hz
        fmin[det.name] = 1.0 / (3600.0 * conf.hours_per_obs * conf.num_obs)
        fknee[det.name] = det.fknee_Hz
        alpha[det.name] = det.alpha
        net[det.name] = det.net

    noise = tt.AnalyticNoise(rate=rate,
                             fmin=fmin,
                             detectors=[x.name for x in conf.detectors],
                             fknee=fknee,
                             alpha=alpha,
                             NET=net)

    # The distributed time stream data
    data = toast.Data(comm)

    intsamp = None
    if len(intervals) == 1:
        intsamp = intervals[0].last - intervals[0].first + 1
    else:
        # include also the samples in the gap
        intsamp = intervals[1].first - intervals[0].first

    intchunk = int(intsamp / conf.chunks_per_obs)
    distint = [] # type: List[int]
    for _ in range(conf.chunks_per_obs - 1):
        distint.append(intchunk)
    distint.append(intsamp - (intchunk * (conf.chunks_per_obs - 1)))

    distsizes = []
    for _ in intervals:
        distsizes.extend(distint)

    totsamples = np.sum(distsizes)

    # create the single TOD for this observation
    tod = tt.TODSatellite(
        mpicomm=comm.comm_group,
        detectors=detquats,
        samples=totsamples,
        firsttime=0.0,
        rate=conf.sample_rate_Hz,
        spinperiod=conf.spin_period_min,
        spinangle=conf.spin_angle_deg,
        precperiod=conf.prec_period_min,
        precangle=conf.prec_angle_deg,
        sizes=distsizes)

    # Create the (single) observation

    ob = {} # type: Dict[str, Any]
    ob['id'] = 0
    ob['name'] = 'mission'
    ob['tod'] = tod
    ob['intervals'] = intervals
    ob['noise'] = noise

    data.obs.append(ob)

    precquat = tt.slew_precession_axis(nsim=tod.local_samples[1],
                                       firstsamp=tod.local_samples[0],
                                       samplerate=conf.sample_rate_Hz,
                                       degday=DEG_ANGLE_PER_DAY)

    # we set the precession axis now, which will trigger calculation
    # of the boresight pointing.
    tod.set_prec_axis(qprec=precquat)

    # simulate noise

    nse = tt.OpSimNoise(out='noise', realization=0)
    nse.exec(data)

    comm.comm_world.barrier()

    # make a Healpix pointing matrix.  By setting purge_pntg=True,
    # we purge the detector quaternion pointing to save memory.
    # If we ever change this pipeline in a way that needs this
    # pointing at a later stage, we need to set this to False
    # and run at higher concurrency.

    pointing = tt.OpPointingHpix(nside=conf.nside,
                                 nest=False,
                                 mode='IQU',
                                 hwprpm=conf.hwp_rpm,
                                 hwpstep=None,
                                 hwpsteptime=None)
    pointing.exec(data)

    comm.comm_world.barrier()

    if comm.comm_world.rank == 0:
        log.info('reading map "%s"', conf.sky_map)
    foreground_map = healpy.read_map(conf.sky_map, 
                                     field=(conf.i_stokes, conf.q_stokes, conf.u_stokes), 
                                     verbose=False)

    for det in conf.detectors:
        if comm.comm_world.rank == 0:
            log.info('processing detector "%s"', det.name)
        times = tod.read_times()
        flags, glflags = tod.read_flags(detector=det.name)

        pixels = tod.cache.reference('pixels_{0}'.format(det.name))
        theta, phi = healpy.pix2ang(conf.nside, pixels)
        vect = healpy.pix2vec(conf.nside, pixels)
        weights = tod.cache.reference('weights_{0}'.format(det.name))
        psi = 0.5 * np.arctan2(weights[:, 2], weights[:, 1])
        psi += np.deg2rad(det.polarisation_angle_deg)

        noise = tod.cache.reference('noise_{0}'.format(det.name))

        # Estimate the total TOD using the formula
        #   P = I + Q cos 2psi + U sin 2psi
        # where Q and U are expressed in the reference frame of the detector
        sky_i, sky_q, sky_u = [healpy.get_interp_val(foreground_map[x], theta, phi)
                               for x in (0, 1, 2)]
        det_q = sky_q * np.cos(2. * psi) + sky_u * np.sin(2. * psi)
        det_u = -sky_q * np.sin(2. * psi) + sky_u * np.cos(2. * psi)
        foreground_tod = sky_i + det_q * np.cos(2. * psi) + det_u * np.sin(2. * psi)

        dipole_tod = dip.get_dipole_temperature(np.column_stack(vect).T)
        total_tod = det.gain * (foreground_tod + dipole_tod + noise)
        hdu = fits.BinTableHDU.from_columns([fits.Column(name='TIME',
                                                         format='D',
                                                         unit='s',
                                                         array=times),
                                             fits.Column(name='THETA',
                                                         format='E',
                                                         unit='rad',
                                                         array=theta),
                                             fits.Column(name='PHI',
                                                         format='E',
                                                         unit='rad',
                                                         array=phi),
                                             fits.Column(name='PSI',
                                                         format='E',
                                                         unit='rad',
                                                         array=psi),
                                             fits.Column(name='FGTOD',
                                                         format='E',
                                                         unit='K',
                                                         array=foreground_tod),
                                             fits.Column(name='DIPTOD',
                                                         format='E',
                                                         unit='K',
                                                         array=dipole_tod),
                                             fits.Column(name='NOISETOD',
                                                         format='E',
                                                         unit='K',
                                                         array=noise),
                                             fits.Column(name='TOTALTOD',
                                                         format='E',
                                                         unit='V',
                                                         array=total_tod),
                                             fits.Column(name='FLAGS',
                                                         format='I',
                                                         unit='',
                                                         array=flags),
                                             fits.Column(name='GLFLAGS',
                                                         format='I',
                                                         unit='',
                                                         array=glflags)])
        hdu.header['COMMENT'] = 'Angles are expressed in the Ecliptic system'
        hdu.header['FIRSTT'] = (times[0], 'Time of the first sample [s]')
        hdu.header['LASTT'] = (times[-1], 'Time of the last sample [s]')
        hdu.header['GAIN'] = (det.gain, 'Detector gain [V/K]')
        hdu.header['NAME'] = (det.name, 'Name of the detector')
        hdu.header['VERSION'] = (__version__,
                                 'Version of the code used to create '
                                 'this file')

        hdu.header['net'] = (det.net, 'net')
        hdu.header['ALPHA'] = (det.alpha, 'Slope of the 1/f noise')
        hdu.header['FKNEE'] = (det.fknee_Hz,
                               'Knee frequency of the 1/f noise [Hz]')

        output_file = os.path.join(outdir,
                                   ('tod_{0}_mpi{1:04d}.fits'
                                    .format(det.name, comm.comm_world.rank)))

        # Save a copy of the parameter file into the primary HDU of the file,
        # so that it will be always possible to know the value of the input parameters
        # used to create it
        prim_hdu = fits.PrimaryHDU(data=conf.param_file_contents)
        hdu_list = fits.HDUList([prim_hdu, hdu])
        hdu_list.writeto(output_file, clobber=True)

    if comm.comm_world.rank == 0:
        log.info('pointings have been created successfully')

##############################################################################

if __name__ == "__main__":
    if len(sys.argv) != 3:
        print(__doc__.format(basename=os.path.basename(sys.argv[0])))
        sys.exit(1)

    main(parameter_file=sys.argv[1], outdir=sys.argv[2])
