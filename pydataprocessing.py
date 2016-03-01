"""This file contains classes and methods for reducing data produced by
PyFIIns.
"""

from __future__ import absolute_import

import astropy.io.fits as pyfits
import collections
import datetime
import os.path
import pp

import common.commonobjects as co
import cleanimage
import dirtyimage
import reduceinterferogram
import renderer
import readfits
import writefits


class PyDataProcessing(object):
    """Class to reduce data from FISICA simulator.
 
    Contains methods:
    __init__
    reduce
    __repr__
    """

    def __init__(self, fitsfile):
        """Constructor.

        Keyword parameters:
        fitsfile     -- name of FITS file with data
        """ 
        self.fitsfile = fitsfile

        self.result = collections.OrderedDict()

        # start parallel python (pp), find the number of CPUS available
        ppservers = ()
        self.job_server = pp.Server(ppservers=ppservers)
        print 'PyFIInS starting pp with %s workers' % \
          self.job_server.get_ncpus()

    def reduce(self, data_quality='noisy', write_idi=True, make_dirty=False):
        """Method invoked to do the work.

        Keyword parameters:
        data_quality -- type of data to be reduced:
                        'pure' will take simulated data with
                          background/detector noise not added
                        'noisy' will take simulated data with added noise
        write_idi    -- True if want to write u-v results to IDI
                        FITS file (for input to CASA).
        make_dirty   -- True if want to form dirty imge from u-v results.

        """
        if data_quality not in ['noisy', 'pure']:
            raise Exception, 'bad data_quality: %s' % data_quality

        # store data and time of run
        now = datetime.datetime.today()
        self.result['runtime'] = now.strftime('%Y%m%dT%H%M%S')

        # read in the interferograms from the FITS file
        fits = readfits.ReadFITS(self.fitsfile)
        self.result['readfits'] = fits.run()   
        print fits

        # recover spectra from interferograms
        reduceint = reduceinterferogram.ReduceInterferogram(
          previous_results=self.result, data_quality=data_quality,
          job_server=self.job_server)
        self.result['reduceinterferogram'] = reduceint.run()
        print reduceint

        if write_idi:
            import writefits_idi

            # write uv data to FITS-IDI
            fitsidi = writefits_idi.WriteFITS_IDI(previous_results=self.result)
            self.result['fits-idi'] = fitsidi.run()
            print fitsidi

        if make_dirty:
            # construct dirty image
            dirty = dirtyimage.DirtyImage(previous_results=self.result,
              job_server=self.job_server)
            self.result['dirtyimage'] = dirty.run()
            print dirty

            # write the cube to FITS
            dirty_image = self.result['dirtyimage']['dirtyimage']
            spatial_axis = self.result['dirtyimage']['spatial axis [arcsec]']
            wn_axis = self.result['dirtyimage']['wavenumber [cm-1]']

            axis1 = co.Axis(data=wn_axis, title='Frequency', units='cm-1')
            axis2 = co.Axis(data=spatial_axis, title='Dec', units='arcsec')
            axis3 = co.Axis(data=spatial_axis, title='RA', units='arcsec')
            cube = co.Image(data=dirty_image, axes=[axis1, axis2, axis3])

            dirtyfitsfile = os.path.splitext(self.fitsfile)[0]
            dirtyfitsfile = '%s.%s.dirtycube.fits' % (dirtyfitsfile,
              data_quality)

            cubewriter = writefits.WriteFITSCube(cube, dirtyfitsfile)
            cubewriter.run()
   
        # construct html description of result
        htmlrenderer = renderer.Renderer(result=self.result)
        htmlrenderer.run(prefix='dp')

    def __repr__(self):
        return 'FISICA Data Processing'

