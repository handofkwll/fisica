from __future__ import absolute_import

import collections
import math
import numpy as np
import astropy.io.fits as pyfits

import common.commonobjects as co


class SkyLoader(object):
    """Class to load a model sky from a FITS file.
    """

    def __init__(self, parameters, previous_results):
        self.parameters = parameters
        self.previous_results = previous_results
        self.result = collections.OrderedDict()

    def run(self):
        print 'SkyLoader.run'

        loadparam = self.previous_results['loadparameters']
        sky_fits = loadparam['substages']['Sky']['FITS file'][0]

        hdulist = pyfits.open(sky_fits)
        print hdulist.info()
        prihdr = hdulist[0].header
        print repr(prihdr)
        skydata = hdulist[0].data

        # looks to be an inversion in index order between FITS and Python
        datashape = np.shape(skydata)
        print datashape
        nx = datashape[2]
        ny = datashape[1]
        if nx != ny:
            raise Exception, 'nx does not equal ny'

        nfreq = datashape[0]

        # pixel increments in degrees
        cdelt1 = prihdr['CDELT1']
        pixsize_rad = np.deg2rad(cdelt1)
        cdelt2 = prihdr['CDELT2']
        cdelt3 = prihdr['CDELT3']

        crval1 = prihdr['CRVAL1']
        crval2 = prihdr['CRVAL2']
        crval3 = prihdr['CRVAL3']

        crpix1 = prihdr['CRPIX1']
        crpix2 = prihdr['CRPIX2']
        crpix3 = prihdr['CRPIX3']

        cunit1 = prihdr['CUNIT1']
        cunit2 = prihdr['CUNIT2']
        cunit3 = prihdr['CUNIT3']

        bunit = prihdr['BUNIT']

        hdulist.close()

        # spatial axes same for all wavelengths - arcsec
        if 'DEG' in cunit1.upper():
            spatial_axis = (crval1 + 
              (np.arange(nx) + 1 - crpix1) * cdelt1) * 3600.0
        else:
            raise Exception, 'cannot handle CUNIT1=%s' % cunit1

        # frequency axis
        if 'HZ' in cunit3.upper():
            # ..in Hz
            fvec = (np.arange(nfreq) - crpix3) * cdelt3 + crval3
            # ..in cm-1
            fvec_wn = fvec / 3.0e10
            # ..wavelength in metres
            fvec_wl = 3.0e8 / fvec
        else:
            raise Exception, 'cannot handle CUNIT3=%s' % cunit1

        # get the frequency range being observed by the FTS
        fts = self.previous_results['fts']
        fts_wn_truncated = fts['fts_wn_truncated']

        skymodel = np.zeros([len(fts_wn_truncated),ny,nx], np.float)

        # interpolate in the spectral dimension from the input model
        for i in range(nx):
            for j in range(ny):
                skymodel[:,j,i] = np.interp(fts_wn_truncated, fvec_wn,
                  skydata[:,j,i], 0.0, 0.0)

        
        # copying Matlab version, remove -ve numbers (noise) as unphysical
        # ..now commented out as +ve noise leads to big I1 + I2 signal 
        # ..which is unrealistic
#        skymodel[skymodel<0.0] = 0.0

        # convert from Jy/pixel to W/m2/cm-1/pixel
        if 'JY/PIXEL' in bunit.upper():
            skymodel *= (1.0e-26 * 299792458.0 * 100) 
        else:
            print 'cannot handle BUNIT=%s' % bunit

#        print 'skymodel has single spike'
#        skymodel[skymodel > 0.0] = 0.0
#        skymodel[:,nx/2,nx/2] = 1.0

        self.result['sky model'] = skymodel
        self.result['spatial axis [arcsec]'] = spatial_axis
        self.result['frequency axis'] = fts_wn_truncated 
        self.result['fvec_wn'] = fvec_wn
        self.result['npix'] = nx
        self.result['pixsize [rad]'] = pixsize_rad

        return self.result

    def __repr__(self):
        return '''
SkyLoader:
  nx, ny, nf  : {shape}
  wnmin       : {wnmin}
  wnmax       : {wnmax}
  native wnmin: {native_wnmin}
  native wnmax: {native_wnmax}
'''.format(
          shape=np.shape(self.result['sky model']),
          wnmin=np.min(self.result['frequency axis']),
          wnmax=np.max(self.result['frequency axis']),
          native_wnmin=np.min(self.result['fvec_wn']),
          native_wnmax=np.max(self.result['fvec_wn']))

