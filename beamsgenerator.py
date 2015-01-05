from __future__ import absolute_import

import collections
import numpy as np

import common.commonobjects as co

def calculate_dirty_beam(npix, rpix, mult, bmax, bxbylist, wnmax, wn):
    lamb = 1.0 / (wn * 100.0)
    lambda_min = 1.0 / (wnmax * 100.0)
    maxbx = bmax / lambda_min

    uv = numpy.zeros([mult*npix,mult*npix])

    for bxby in bxbylist:
        ibx = int(round(((bxby[0] / lamb) / maxbx) * mult * rpix))
        if ibx < 0:
            ibx = mult*npix + ibx 
        iby = int(round(((bxby[1] / lamb) / maxbx) * mult * rpix))
        if iby < 0:
            iby = mult*npix + iby 
        uv[ibx,iby] = 1.0
 
    temp = numpy.fft.fftshift(uv)
    temp = numpy.fft.fft2(temp)
    temp = numpy.fft.fftshift(temp)
    temp = numpy.power(temp, 2)
    # convert to real, should be real anway as FT of symmetric function 
    dirty_beam = numpy.abs(temp)
    # and normalise
    dirty_beam /= numpy.max(dirty_beam)

    # truncate dirty beam pattern
    dirty_beam = dirty_beam[
      (mult-1)*rpix:(mult+1)*rpix, (mult-1)*rpix:(mult+1)*rpix]

    # fit Gaussian to dirty beam
    shape = numpy.shape(dirty_beam)
    fitter = fitgaussian.FitGaussian()
    dirty_centre = numpy.array(dirty_beam[shape[0]/2-5:shape[0]/2+5,
      shape[1]/2-5:shape[1]/2+5])
    p = fitter.fitgaussian(dirty_centre)

    # clean beam
    cp = (1.0, float(shape[0])/2.0, float(shape[1])/2.0, p[3], p[4], p[5])
    rotgauss = fitter.gaussian(*cp)
    bmax = numpy.max(dirty_beam)
    clean_beam = numpy.fromfunction(rotgauss, numpy.shape(dirty_beam))

    return dirty_beam, clean_beam

def calculate_primary_beam(npix, rpix, mult, mx, my, bmax, m1_diameter,
  wnmax, wn):
    # assume uniform illumination of primary. uv goes out to
    # +- bmax / lambda_min, fill uv plane appropriate to primary
    # size
    mirror = numpy.ones([mult*npix,mult*npix])
    lambda_min = 1.0 / (wnmax * 100.0)
    lamb = 1.0 / (wn * 100.0)
    mirror[(numpy.sqrt(mx*mx + my*my) / (mult*rpix)) * 
      (bmax / lambda_min) > (0.5 * m1_diameter / lamb)] = 0.0

    # fftshift mirror to be centred at [0,0],
    # do a 2d ftt of it,
    # fftshift transform so that 0 freq is in centre of array
    temp = numpy.fft.fftshift(mirror)
    temp = numpy.fft.fft2(temp)
    temp = numpy.fft.fftshift(temp)
    # convert to real, should be real anyway as FT of symmetric function 
    # square to give intensity response
    primary_amplitude_beam = numpy.real(temp)
    primary_intensity_beam = numpy.power(numpy.abs(temp), 2)
    # and normalise
    primary_amplitude_beam /= numpy.max(primary_amplitude_beam)
    primary_intensity_beam /= numpy.max(primary_intensity_beam)

    # truncate primary beam pattern
    primary_amplitude_beam = primary_amplitude_beam[
      (mult-1)*rpix:(mult+1)*rpix,
      (mult-1)*rpix:(mult+1)*rpix]
    primary_intensity_beam = primary_intensity_beam[
      (mult-1)*rpix:(mult+1)*rpix,
      (mult-1)*rpix:(mult+1)*rpix]

    return primary_intensity_beam, primary_amplitude_beam


class BeamsGenerator(object):
    """Class to generate the primary beam(s) of the simulated observation.
    """

    def __init__(self, previous_results, job_server):
        self.previous_results = previous_results
        self.result = collections.OrderedDict()
        self.job_server = job_server

    def run(self):
        print 'Calculating beams'

        # gather configuration 
        uvmapgen = self.previous_results['uvmapgenerator']
        fts = self.previous_results['fts']
        telescope = self.previous_results['loadparameters']['substages']\
          ['Telescope']

        # number of mirror positions sampled by FTS
        nsample = fts['ftsnsample']
        fts_wn = fts['fts_wn']
        fts_wn_truncated = fts['fts_wn_truncated']
      
        # max pixel size (radians) that will fully sample the image.
        # Sampling freq is 2 * Nyquist freq.
        # So sampling freq = 2 * b_max / lambda_min.
        bmax = uvmapgen['bmax']
        self.result['max baseline [m]'] = bmax

        lambda_min = 1.0 / (np.max(fts_wn) * 100.0)
        max_pixsize = lambda_min / (2.0 * bmax)
        self.result['max pixsize [rad]'] = max_pixsize

        # Number of pixels per beam
        # radius of first null of Airy disk is theta=1.22lambda/d. Number of
        # pixels is beam diameter/max_pixsize. Calculate this for the 
        # longest wavelength that has flux (at wnmin), largest beam
        
        max_beam_radius = 1.22 * (1.0 / (np.min(fts_wn_truncated) * 100.0)) /\
          telescope['Primary mirror diameter']
        self.result['beam diam [rad]'] = 2.0 * max_beam_radius
         
        # go out to radius of first null
        rpix = int(max_beam_radius / max_pixsize)
        npix = 2 * rpix
        self.result['npix'] = npix

        # spatial axes same for all wavelengths
        axis = np.arange(-rpix, rpix, dtype=np.float)
        axis *= max_pixsize
        axis *= 206264.806247
        self.result['spatial axis [arcsec]'] = axis
        axis1 = co.Axis(data=-axis, title='RA offset', units='arcsec') 
        axis2 = co.Axis(data=axis, title='Dec offset', units='arcsec') 
        axis3 = co.Axis(data=fts_wn_truncated, title='Frequency', units='cm-1')

        # calculate beams for each point on fts_wn_truncated
        # oversample uv grid by mult to minimise aliasing
        mult = 5
        mx,my = np.meshgrid(np.arange(-mult*rpix, mult*rpix),
          np.arange(-mult*rpix, mult*rpix))
        self.result['primary beam'] = collections.OrderedDict()
        self.result['primary amplitude beam'] = collections.OrderedDict()
        self.result['dirty beam'] = collections.OrderedDict()
        self.result['clean beam'] = collections.OrderedDict()

        wnmax = np.max(fts_wn_truncated)
        m1_diameter = telescope['Primary mirror diameter']

        # first calculate primary beam for each wn
        jobs = {}
        for wn in fts_wn_truncated:
            # submit jobs
            indata = (npix, rpix, mult, mx, my, bmax, m1_diameter, wnmax,
              wn,)
            jobs[wn] = self.job_server.submit(calculate_primary_beam,
              indata, (), ('numpy',))

        primary_amplitude_beam = np.zeros([npix,npix,len(fts_wn_truncated)],
          np.complex)
        for iwn,wn in enumerate(fts_wn_truncated):
            # collect and store results
            primary_intensity_beam, primary_amplitude_beam[:,:,iwn] = jobs[wn]()
            image = co.Image(data=primary_intensity_beam, axes=[axis1, axis2],
              title='Primary Beam %06.4g cm-1' % wn)
            self.result['primary beam'][wn] = image

        cube = co.Cube(data=primary_amplitude_beam, axes=[axis1, axis2, axis3],
          title='Amplitude Primary Beam %06.4g cm-1' % wn)
        self.result['primary amplitude beam'] = cube

        # second calculate dirty beam, same use of pp
        for wn in fts_wn_truncated:
            #
            # 1. uv plane covered by discrete points so should use a
            # slow dft rather than fft, which requires evenly spaced
            # samples. Perhaps a bit slow though as calcs done in Python
            # layer.
            # 2. Could use shift theorem to give even-spaced equivalent
            # of each discrete point, then use ffts. Equivalent to 1
            # in speed.
            # 3. Oversample uv plane and accept small error in point
            # positions. Fastest but adequate?

            indata = (npix, rpix, mult, bmax, uvmapgen['bxby'], wnmax, wn,)
            jobs[wn] = self.job_server.submit(calculate_dirty_beam,
              indata, (), ('numpy','fitgaussian',))


        # axes for dirty beam
        axis = np.arange(-rpix, rpix, dtype=np.float)
        axis *= (lambda_min / bmax)
        axis *= 206264.806247
        axis1 = co.Axis(data=-axis, title='RA offset', units='arcsec') 
        axis2 = co.Axis(data=axis, title='Dec offset', units='arcsec') 

        for wn in fts_wn_truncated:
            dirty_beam,clean_beam = jobs[wn]()
            image = co.Image(data=dirty_beam, axes=[axis1, axis2],
              title='Dirty Beam %06.4g cm-1' % wn)
            self.result['dirty beam'][wn] = image
            image = co.Image(data=clean_beam, axes=[axis1, axis2],
              title='Clean Beam %06.4g cm-1' % wn)
            self.result['clean beam'][wn] = image

        return self.result

    def __repr__(self):
        return 'PrimaryBeamGenerator'

