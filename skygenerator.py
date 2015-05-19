from __future__ import absolute_import

import collections
import math
import numpy as np
import scipy.constants as sc

import common.commonobjects as co

def black_body(temperature, wavenumber):
    """Function to calculate Planck function. 
       temperature - Kelvin
       wavenumber - frequency in cm-1
    """
    freq = wavenumber * sc.c * 100.0
    if freq > 0:
        jnu = 2.0 * sc.h * pow(freq,3) / (pow(sc.c,2) *
          (math.exp((sc.h * freq) / (sc.k * temperature)) - 1.0))
    else:
        jnu = 0.0

    return jnu


class BB_spectrum(object):
    """Class to generate BB spectrum.
    """
    def __init__(self, temperature, frequency_axis, cutoffmin, cutoffmax,
      emissivity):
        self.temperature = temperature
        self.frequency_axis = frequency_axis
        self.cutoffmin = cutoffmin
        self.cutoffmax = cutoffmax
        self.emissivity = emissivity

    def calculate(self):
        spectrum = np.zeros(np.shape(self.frequency_axis))
        # ignore floating point errors
        old_settings = np.seterr(all='ignore')

        for iwn,wn in enumerate(self.frequency_axis):
            # simulate real-life 'rounded' cutoffs numerically  
            f1 = 1.0 / (1.0 + pow(self.cutoffmin/wn, 18) + 
              pow(wn/self.cutoffmax, 24))
            spectrum[iwn] = black_body(self.temperature, wn) * f1 * \
              self.emissivity

        # make sure spectrum is zero beyond cutoffs, can cause 
        # aliasing at high freq
        spectrum[self.frequency_axis <= self.cutoffmin] = 0.0
        spectrum[self.frequency_axis >= self.cutoffmax] = 0.0

        # restore fp behaviour
        ignore = np.seterr(**old_settings)

        return spectrum


class linespectrum(object):
    """Class to return spectrum with line at wn.
    """
    def __init__(self, temperature, linefreq, frequency_axis):
        self.temperature = temperature
        self.linefreq = linefreq
        self.frequency_axis = frequency_axis

    def calculate(self):
        # spectrum is a sync function centred on the specified line freq
        x = (self.frequency_axis - self.linefreq) / \
          (self.frequency_axis[1] - self.frequency_axis[0])
        spectrum = np.sinc(x) * self.temperature

        return spectrum


class SkyGenerator(object):
    """Class to generate a model sky.
    """

    def __init__(self, parameters, previous_results):
        self.parameters = parameters
        self.previous_results = previous_results
        self.result = collections.OrderedDict()

    def cubeparameters(self, fts_wn, bmax, m1_diam):
        """Routine to calculate the parameters of the cube holding
        the sky model.
        """
        # max pixel size (radians) that will fully sample the image.
        # Sampling freq is 2 * Nyquist freq.
        # So sampling freq = 2 * b_max / lambda_min.
        lambda_min = 1.0 / (np.max(fts_wn) * 100.0)
        max_pixsize = lambda_min / (2.0 * bmax)
        oversampling = 3.0
        pixsize = max_pixsize / oversampling
#        self.result['pixsize [rad]'] = pixsize = max_pixsize / oversampling
#        self.result['pixsize [arcsec]'] = np.rad2deg(pixsize) * 3600.0

        # Number of pixels per primary beam - measured out to first
        # null of Airy disk.
        # Radius of first null of Airy disk is theta=1.22lambda/d. Number of
        # pixels is beam diameter/max_pixsize. Calculate this for the 
        # longest wavelength that has flux (at wnmin), largest beam
        max_beam_radius = 1.22 * (1.0 / (np.min(fts_wn) * 100.0)) /\
          m1_diam
#        self.result['beam diam [rad]'] = 2.0 * max_beam_radius
         
        # go out to radius of first null
        rpix = int(max_beam_radius / pixsize)
        print 'cubeparameters out of half radius of first null'
        rpix /= 2
        npix = 2 * rpix
#        self.result['npix'] = npix

        # spatial axes same for all wavelengths
        spatial_axis = np.arange(-rpix, rpix, dtype=np.float)
        spatial_axis *= pixsize
        spatial_axis = np.rad2deg(spatial_axis) * 3600.0
#        self.result['spatial axis [arcsec]'] = axis

        return npix, pixsize, spatial_axis

    def run(self):
        print 'SkyGenerator.run'

        # gather configuration info required
        fts = self.previous_results['fts']
        fts_wn_truncated = fts['fts_wn_truncated']
        # truncate spectrum inside allowed spectral range to prevent
        # aliasing problem at wnmax
        cutoffmin = fts['wnmin']
        cutoffmax = fts['wnmax'] - 1.0

        uvmapgen = self.previous_results['uvmapgenerator']
        bmax = uvmapgen['bmax']

        telescope = self.previous_results['loadparameters']['substages']\
          ['Telescope']
        m1_diam = telescope['Primary mirror diameter']

        # calculate the spatial dimensions of the cube
        npix, pixsize, spatial_axis = self.cubeparameters(fts_wn_truncated,
          bmax, m1_diam)

        # skymodel is complex so that its fft can hold truncated version
        # of infinitesimally sampled map - does that make sense?
        skymodel = np.zeros([npix, npix, len(fts_wn_truncated)], np.complex)

        sky = self.parameters['substages']['Sky']
        columns = sky['sourcenum'].keys()

        self.result['sources'] = collections.OrderedDict()

        for column in columns:
            temp = sky['sourcenum'][column]
            sourcenum = int(round(temp))
            if sourcenum not in self.result['sources'].keys():
                self.result['sources'][sourcenum] = {}

            type = sky['type'][column]

            temp = sky['x pos [asec]'][column]
            xpos = float(temp)

            temp = sky['y pos [asec]'][column]
            ypos = float(temp)

            temp = sky['xwidth'][column]
            try:
                xwidth = float(temp)
            except:
                xwidth = None

            temp = sky['ywidth'][column]
            try:
                ywidth = float(temp)
            except:
                ywidth = None

            spectrum = sky['spectrum'][column]

            temp = sky['temperature'][column]
            temperature = float(temp)

            temp = sky['linefreq'][column]
            try:
                linefreq = float(temp)
            except:
                linefreq = None

            temp = sky['emissivity'][column]
            emissivity = float(temp)

            print 'generating source:%s type:%s xpos:%s ypos:%s' % (sourcenum, 
              type, xpos, ypos)
            if type.upper().strip() == 'GAUSSIAN':
                print '  fwhmx:%s fwhmy:%s' % (xwidth, ywidth)

            if spectrum.upper().strip() == 'BLACKBODY':
                print '  blackbody spectrum temperature:%s cutoffmin:%s cutoffmax:%s e:%s' % (
                  temperature, cutoffmin, cutoffmax, emissivity)
                spectrum_func = BB_spectrum(temperature, fts_wn_truncated,
                  cutoffmin, cutoffmax, emissivity)

            elif spectrum.upper().strip() == 'LINE':
                print '  line spectrum brightness temperature:%s cutoffmin:%s cutoffmax:%s e:%s' % (
                  temperature, cutoffmin, cutoffmax, emissivity)
                spectrum_func = linespectrum(temperature, linefreq, fts_wn_truncated)

            if type.upper().strip() == 'POINT':
                source_spectrum = self._create_point_source(xpos, ypos,
                  skymodel, spatial_axis, fts_wn_truncated, spectrum_func)
            elif type.upper().strip() == 'GAUSSIAN':
                source_spectrum = self._create_gaussian_source(xpos, ypos,
                  xwidth, ywidth, skymodel, spatial_axis, fts_wn_truncated,
                  spectrum_func)
            else:
                source_spectrum = None
                print "source type '%s' not yet implemented" % type

            self.result['sources'][sourcenum]['spectrum'] = source_spectrum

        self.result['sky model'] = skymodel
        self.result['spatial axis [arcsec]'] = spatial_axis
        self.result['frequency axis'] = fts_wn_truncated 
        self.result['pixsize [rad]'] = pixsize

        return self.result

    def _create_gaussian_source(self, xpos, ypos, xwidth, ywidth,
      skymodel, spatial_axis, frequency_axis, spectrum_function):
        """Create a point source.
        """
        # calculate spectrum
        spectrum = spectrum_function.calculate()

        # calculate Gaussian profile - a cunning use of array slicing found
        # on the web
        x = spatial_axis
        y = x[:,np.newaxis]
        profile = np.exp(-4*np.log(2) * ((x-xpos)**2/xwidth**2 + (y-ypos)**2/ywidth**2))

        # apodise
        profile[((x-xpos)**2 + (y-ypos)**2) > 16] = 0.0

        # go through freq planes and add source to each
        for iwn,wn in enumerate(frequency_axis):
            # add to sky model
            skymodel[:,:,iwn] += profile * spectrum[iwn]

        # return spectrum
        axis = co.Axis(data=frequency_axis, title='wavenumber', units='cm-1')
        spectrum = co.Spectrum(data=spectrum, axis=axis,
          title='Source spectrum', units='W sr-1 m-2 Hz-1')

        return spectrum

    def _create_point_source(self, xpos, ypos, skymodel, spatial_axis,
      frequency_axis, spectrum_function):
        """Create a point source.
        """

        # calculate xpos, ypos in units of pixel - numpy arrays [row,col]
        nx = len(spatial_axis)
        colpos = float(nx-1) * float (xpos - spatial_axis[0]) / (spatial_axis[-1] - spatial_axis[0])
        rowpos = float(nx-1) * float (ypos - spatial_axis[0]) / (spatial_axis[-1] - spatial_axis[0])

        if colpos < 0 or colpos > (nx-1) or rowpos < 0 or rowpos > (nx-1):
            # point source is outside modelled area 
            return

        # calculate fourier phase shift to move point at [0,0] to 
        # [rowpos, colpos]
        shiftx = np.zeros([nx], np.complex)
        shiftx[:nx/2] = np.arange(nx/2, dtype=np.complex)
        shiftx[nx/2:] = np.arange(-nx/2, 0, dtype=np.complex)
        shiftx = np.exp((-2.0j * np.pi * colpos * shiftx) / float(nx))

        shifty = np.zeros([nx], np.complex)
        shifty[:nx/2] = np.arange(nx/2, dtype=np.complex)
        shifty[nx/2:] = np.arange(-nx/2, 0, dtype=np.complex)
        shifty = np.exp((-2.0j * np.pi * rowpos * shifty) / float(nx))

        shift = np.ones([nx,nx], np.complex)
        for j in range(nx):
            shift[j,:] *= shiftx
        for i in range(nx):
            shift[:,i] *= shifty

        # calculate spectrum
        spectrum = spectrum_function.calculate()

        # go through freq planes and add point source to each
        for iwn,wn in enumerate(frequency_axis):
            # create point in frequency space
            temp = np.zeros([nx,nx])
            temp[0,0] = spectrum[iwn]
            # 2d fft
            temp = np.fft.fft2(temp)
            # apply phase shift to move point to required offset
            temp *= shift
            # transform back
            temp = np.fft.ifft2(temp)

            # add to sky model
            skymodel[:,:,iwn] += temp

        # return spectrum
        axis = co.Axis(data=frequency_axis, title='wavenumber', units='cm-1')
        spectrum = co.Spectrum(data=spectrum, axis=axis,
          title='Source spectrum', units='W sr-1 m-2 Hz-1')

        return spectrum

    def __repr__(self):
        return '''
SkyGenerator:
  wnmin    : {wnmin} [cm-1]
  wnmax    : {wnmax} [cm-1]
  delta wn : {delta_wn} [cm-1]
  num wn   : {nwn}
  pixsize  : {pixsize} [arcsec]
  npix     : {npix}
'''.format(
          wnmin = min(self.result['frequency axis']),
          wnmax = max(self.result['frequency axis']),
          delta_wn = abs(self.result['frequency axis'][1] -
          self.result['frequency axis'][0]),
          nwn = len(self.result['frequency axis']),
          pixsize = abs(self.result['spatial axis [arcsec]'][1] - 
          self.result['spatial axis [arcsec]'][0]),
          npix = len(self.result['spatial axis [arcsec]']))

