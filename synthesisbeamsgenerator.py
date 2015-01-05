from __future__ import absolute_import

import collections
import numpy as np

import common.commonobjects as co

def calculate_dirty_beam(npix, pixsize, wn, mult, bxbylist):
    lamb = 1.0 / (wn * 100.0)

    # maximum baseline in the UV plane corresponds to the cube pixsize
    # at this wavelength
    maxbx = lamb / (2.0 * pixsize)

    # oversample UV plane by factor 'mult' to reduce aliasing    
    uv = numpy.zeros([mult*npix, mult*npix])

    # populate uv plane with measured baselines
    rpix = npix / 2
    for bxby in bxbylist:
        ibx = int(round((bxby[0] / maxbx) * mult * rpix))
        if ibx < 0:
            ibx = mult*npix + ibx 
        iby = int(round((bxby[1] / maxbx) * mult * rpix))
        if iby < 0:
            iby = mult*npix + iby 
        uv[ibx,iby] = 1.0 

    # fft the uv plane to get the dirty beam, shift beam centre from
    # 0,0 to array centre
    temp = numpy.fft.fft2(uv)
    temp = numpy.fft.fftshift(temp)

    # calculate normalised intensity beam
    temp *= numpy.conjugate(temp)
    dirty_beam = numpy.abs(temp)
    dirty_beam /= numpy.max(dirty_beam)

    # truncate beam pattern
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


class SynthesisBeamsGenerator(object):
    """Class to generate the dirty and clean beams of the simulated observation.
    """

    def __init__(self, previous_results, job_server):
        self.previous_results = previous_results
        self.result = collections.OrderedDict()
        self.job_server = job_server

    def run(self):
        print 'Calculating synthesis dirty beam...'

        # gather configuration 
        uvmapgen = self.previous_results['uvmapgenerator']
        cubeparams = self.previous_results['cubeparameters']
        self.result['wn'] = wn = cubeparams['wn']
        self.result['pixsize [rad]'] = pixsize = cubeparams['pixsize [rad]']
        self.result['npix'] = npix = cubeparams['npix']

        # spatial axes same for all wavelengths
        self.result['spatial axis [arcsec]'] = axis = \
          cubeparams['spatial axis [arcsec]']
        axis1 = co.Axis(data=-axis, title='RA offset', units='arcsec') 
        axis2 = co.Axis(data=axis, title='Dec offset', units='arcsec') 
        axis3 = co.Axis(data=wn, title='Frequency', units='cm-1')

        # calculate beam for each wavenumber in 'wn'
        # oversample uv grid by mult to minimise aliasing
        rpix = npix / 2
        mult = 5
        mx,my = np.meshgrid(np.arange(-mult*rpix, mult*rpix),
          np.arange(-mult*rpix, mult*rpix))
        self.result['dirty beam'] = collections.OrderedDict()
        self.result['clean beam'] = collections.OrderedDict()

        # calculate synthesis beams using pp
        jobs = {}
        wnmax = np.max(wn)
        for wavenum in wn:
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
            indata = (npix, pixsize, wavenum, mult, uvmapgen['bxby'],)
            jobs[wavenum] = self.job_server.submit(calculate_dirty_beam,
              indata, (), ('numpy','fitgaussian',))

        for wavenum in wn:
            dirty_beam,clean_beam = jobs[wavenum]()
            image = co.Image(data=dirty_beam, axes=[axis1, axis2],
              title='Dirty Beam %06.4g cm-1' % wavenum)
            self.result['dirty beam'][wavenum] = image
            image = co.Image(data=clean_beam, axes=[axis1, axis2],
              title='Clean Beam %06.4g cm-1' % wavenum)
            self.result['clean beam'][wavenum] = image

        return self.result

    def __repr__(self):
        return 'PrimaryBeamGenerator'

