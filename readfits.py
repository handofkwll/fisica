"""This module contains classes and methods used to read data from a FITS
file holding an observation simulated by PyFIInS.
"""

from __future__ import absolute_import

import collections
import numpy as np
import os.path
import astropy.io.fits as pyfits

import timelinegenerator


class ReadFITS(object):
    """Class to read in an observed timeline from a FITS Table.

    Contains methods:
    __init__
    run
    __repr__
    """

    def __init__(self, fitsfile):
        """Constructor.

        Keyword parameters:
        fitsfile - Name of the FITS file containing the data.
        """
        self.fitsfile = fitsfile
        self.result = {}

    def run(self):
        """Method to read the FITS data.
        """
        #print 'ReadFITS.run'

        hdulist = pyfits.open(self.fitsfile)
        hdulist.info()
        prihdr = hdulist[0].header
        print '..header information'
        print repr(prihdr)

        # general info needed by the reduction
        smec_opd_to_mpd = prihdr['smec_o2m']
        wnmin = prihdr['wnmin']
        spatial_axis = hdulist[1].data

        table = hdulist[2].data
        cols = hdulist[2].columns
        cols.info()

        obs_timeline = {}

        times = table.field('Time')
        scans = table.field('Scan')
        bxs = table.field('Baseline x')
        bys = table.field('Baseline y')
        bzs = table.field('Baseline z')
        smec_positions = table.field('SMEC Position')
        smec_nominal_positions = table.field('SMEC Nominal Position')
        smec_velocity_errors = table.field('SMEC Velocity Error')
        flag = table.field('Flag')
        data = table.field('Data')
        pure_data = table.field('Pure Data')
        cr_data = table.field('Cosmic Ray Data')
        detector_noise_data = table.field('Detector Noise Data')
        pointing1_x = table.field('Pointing1 X')
        pointing1_y = table.field('Pointing1 Y')
        pointing2_x = table.field('Pointing2 X')
        pointing2_y = table.field('Pointing2 Y')
        
        for i,t in enumerate(times):
            obs_timeline[t] = timelinegenerator.Config(
              scans[i], t,
              bxs[i], bys[i], bzs[i], 0,
              smec_positions[i],
              smec_nominal_positions[i],
              flag[i],
              smec_velocity_errors[i],
              pointing1_x[i], pointing1_y[i],
              pointing2_x[i], pointing2_y[i],
              data[i], pure_data[i], cr_data[i], detector_noise_data[i])

        self.result['observed_timeline'] = obs_timeline
        self.result['smec_opd_to_mpd'] = smec_opd_to_mpd
        self.result['wnmin'] = wnmin
        self.result['spatial axis [arcsec]'] = spatial_axis
        self.result['fitsfile'] = self.fitsfile
        return self.result

    def __repr__(self):
        return '''
ReadFITS:
  FITS file : {fitsfile}
'''.format(fitsfile=self.result['fitsfile'])
