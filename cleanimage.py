from __future__ import absolute_import

import collections
import math
import numpy as np
import pp

import common.commonobjects as co
import pythonclean


class CleanImage(object):
    """Class to compute clean image cube from uvspectra.
    """

    def __init__(self, previous_results, job_server):
        self.previous_results = previous_results
        self.job_server = job_server

        self.result = collections.OrderedDict()      

    def run(self):
        print 'CleanImage.run'

        # get dirty image and beam info
        dirty = self.previous_results['dirtyimage']
        dirtyimage = dirty['dirtyimage']
        dirtybeam = dirty['dirtybeam']
        spatial_axis = dirty['spatial axis [arcsec]']
        wavenumber = dirty['wavenumber [cm-1]']

        # clean image cube
        cleanimage = np.zeros(np.shape(dirtyimage), np.float)
        residualimage = np.zeros(np.shape(dirtyimage), np.float)

        # calculate clean image for each wn
        jobs = {}
        # clean only central quarter
        imsize = np.shape(dirtyimage)[:2]
        window = [slice(imsize[0]/4, imsize[0]*3/4), slice(imsize[1]/4, imsize[1]*3/4)]
        for iwn,wn in enumerate(wavenumber):
            # submit jobs
            immax = np.max(dirtyimage[:,:,iwn])
            indata = (dirtyimage[:,:,iwn], dirtybeam[:,:,iwn], window, 0.1,
              immax/100.0, 5000,)
            jobs[wn] = self.job_server.submit(pythonclean.hogbom,
              indata, (), ('numpy','pythonclean',))

        for iwn,wn in enumerate(wavenumber[:3]):
            # collect and store results
            cleanimage[:,:,iwn], residualimage[:,:,iwn] = jobs[wn]()
            print wn, 'clean finished'

        self.result['dirtyimage'] = dirtyimage
        self.result['cleanimage'] = cleanimage + residualimage
        self.result['residualimage'] = residualimage
        self.result['spatial axis [arcsec]'] = spatial_axis
        self.result['spatial axis'] = spatial_axis
        self.result['wavenumber [cm-1]'] = wavenumber

        return self.result
                
    def __repr__(self):
        return '''
CleanImage:
'''.format(
          num_uvspectra=len(self.uvspectra))

