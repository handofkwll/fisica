from __future__ import absolute_import

import astropy.io.fits as pyfits
import collections
import datetime
import pp

import cleanimage
import dirtyimage
import reduceinterferogram
import renderer
import readfits


class PyDataProcessing(object):
    """Data processing module to reduce data from FISICA simulator.
    """
    def __init__(self, fitsfile):
        self.fitsfile = fitsfile

        self.result = collections.OrderedDict()

        # start parallel python (pp), find the number of CPUS available
        ppservers = ()
        self.job_server = pp.Server(ppservers=ppservers)
        print 'PyFIInS starting pp with %s workers' % \
          self.job_server.get_ncpus()

    def reduce(self):
        # store data and time of run
        now = datetime.datetime.today()
        self.result['runtime'] = now.strftime('%Y%m%dT%H%M%S')

        # read in the interferograms from the FITS file
        fits = readfits.ReadFITS(self.fitsfile)
        self.result['readfits'] = fits.run()   
        print fits

        # recover spectra from interferograms
        reduceint = reduceinterferogram.ReduceInterferogram(
          previous_results=self.result,
          job_server=self.job_server)
        self.result['reduceinterferogram'] = reduceint.run()
        print reduceint
   
        # construct dirty image
        dirty = dirtyimage.DirtyImage(
          previous_results=self.result,
          job_server=self.job_server)
        self.result['dirtyimage'] = dirty.run()
        print dirty
   
        # construct clean image
#    #    clean = cleanimage.CleanImage(
#    #      previous_results=self.result,
#    #      job_server=self.job_server)
#    #    self.result['cleanimage'] = clean.run()
#    #    print clean
#    #

        # construct html description of result
        htmlrenderer = renderer.Renderer(result=self.result)
        htmlrenderer.run(prefix='dp')

    def __repr__(self):
        return 'FISICA Data Processing'

