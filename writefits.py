from __future__ import absolute_import

import collections
import numpy as np
import os.path
import astropy.io.fits as pyfits

import matplotlib.pyplot as plt


class WriteFITS(object):
    """Class to write out observed timeline in a FITS Table.
    """

    def __init__(self, previous_results):
        """Constructor.

        Parameters:
        previous_results - simulator result structure.
        """
        self.previous_results = previous_results
        self.result = {}

    def run(self):
        """Method to write the FITS data.
        """
        #print 'WriteFITS.run'

        # construct the name of the file
        runtime = self.previous_results['runtime']
        fitsname = '%s.fits' % runtime

        # get list of instrument observations
        observe = self.previous_results['observe']
        obs_timeline = observe['observed_timeline']
        observed_times = obs_timeline.keys()
        observed_times.sort()

        # construct lists of the values to be stored in each Table column
        for t in observed_times:
            timelist = []
            smec_position = []
            smec_nominal_position = []
            flag = []
            data = []
            pointing1_x = []
            pointing1_y = []
            pointing2_x = []
            pointing2_y = []

            config = obs_timeline[t]

            timelist.append(config.time)
            smec_position.append(config.smec_position)
            smec_nominal_position.append(config.smec_nominal_position)
            flag.append(config.flag)
            data.append(config.data)
            pointing1_x.append(config.pointing1_x)
            pointing1_y.append(config.pointing1_y)
            pointing2_x.append(config.pointing2_x)
            pointing2_y.append(config.pointing2_y)

        # create a Header object and primary HDU - this just contains
        # some very basic, general information
        prihdr = pyfits.Header()
        prihdr['COMMENT'] = 'This FITS file was created by pyfiins at %s' % \
          runtime
        prihdu = pyfits.PrimaryHDU(header=prihdr)

        # create list of Header Data Unit objects, include the primary HDU
        hdulist = pyfits.HDUList([prihdu])

        # create an HDU to contain the Table and append it to the list
        hdulist.append(pyfits.BinTableHDU.from_columns(
          pyfits.ColDefs([
          pyfits.Column(name='Time', format='D',
          array=np.array(timelist)),
          pyfits.Column(name='SMEC Position', format='E',
          array=np.array(smec_position)),
          pyfits.Column(name='SMEC Nominal Position', format='E',
          array=np.array(smec_nominal_position)),
          pyfits.Column(name='Flag', format='L',
          array=np.array(flag)),
          pyfits.Column(name='Data', format='E',
          array=np.array(data)),
          pyfits.Column(name='Pointing1 X', format='E',
          array=np.array(pointing1_x)),
          pyfits.Column(name='Pointing1 Y', format='E',
          array=np.array(pointing1_y)),
          pyfits.Column(name='Pointing2 X', format='E',
          array=np.array(pointing2_x)),
          pyfits.Column(name='Pointing2 Y', format='E',
          array=np.array(pointing2_y))])))

        # write the HDU list to a file
        hdulist.writeto(fitsname, clobber=True)
        self.result['fitsfile'] = fitsname

        return self.result

    def __repr__(self):
        return '''
WriteFITS:
  FITS file : {fitsfile}
'''.format(fitsfile=self.result['fitsfile'])
