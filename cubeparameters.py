from __future__ import absolute_import

import collections
import numpy as np


class CubeParameters(object):
    """Class to generate the parameters of the simulated cube.
    """

    def __init__(self, previous_results):
        self.previous_results = previous_results
        self.result = collections.OrderedDict()

    def run(self):
        #print 'Calculating cube parameters...'

        # gather configuration info required
        uvmapgen = self.previous_results['uvmapgenerator']
        fts = self.previous_results['fts']
        telescope = self.previous_results['loadparameters']['substages']\
          ['Telescope']

        # wavenumbers of spectrum corresponding to FTS configuration
        self.result['wn'] = fts_wn = fts['fts_wn_truncated']
      
        # max pixel size (radians) that will fully sample the image.
        # Sampling freq is 2 * Nyquist freq.
        # So sampling freq = 2 * b_max / lambda_min.
        bmax = uvmapgen['bmax']
        lambda_min = 1.0 / (np.max(fts_wn) * 100.0)
        max_pixsize = lambda_min / (2.0 * bmax)
        oversampling = 3.0
        self.result['pixsize [rad]'] = pixsize = max_pixsize / oversampling
        self.result['pixsize [arcsec]'] = np.rad2deg(pixsize) * 3600.0

        # Number of pixels per primary beam - measured out to first
        # null of Airy disk.
        # Radius of first null of Airy disk is theta=1.22lambda/d. Number of
        # pixels is beam diameter/max_pixsize. Calculate this for the 
        # longest wavelength that has flux (at wnmin), largest beam
        max_beam_radius = 1.22 * (1.0 / (np.min(fts_wn) * 100.0)) /\
          telescope['Primary mirror diameter']
        self.result['beam diam [rad]'] = 2.0 * max_beam_radius
         
        # go out to radius of first null
        rpix = int(max_beam_radius / pixsize)
        npix = 2 * rpix
        self.result['npix'] = npix

        # spatial axes same for all wavelengths
        axis = np.arange(-rpix, rpix, dtype=np.float)
        axis *= pixsize
        axis = np.rad2deg(axis) * 3600.0
        self.result['spatial axis [arcsec]'] = axis

        return self.result

    def __repr__(self):
        return '''
CubeParameters:
  wnmin    : {wnmin} [cm-1]
  wnmax    : {wnmax} [cm-1]
  delta wn : {delta_wn} [cm-1]
  num wn   : {nwn}
  pixsize  : {pixsize} [arcsec]
  npix     : {npix}
'''.format(
          wnmin = min(self.result['wn']),
          wnmax = max(self.result['wn']),
          delta_wn = abs(self.result['wn'][1] - self.result['wn'][0]),
          nwn = len(self.result['wn']),
          pixsize = self.result['pixsize [arcsec]'],
          npix = self.result['npix'])
