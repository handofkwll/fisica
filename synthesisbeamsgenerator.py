from __future__ import absolute_import

import collections
import numpy as np

import common.commonobjects as co

def calculate_synthesis_beams(npix, pixsize, wn, bxbylist):
    """Routine to calculate the dirty and clean beams on a 
    sky map with npix x npix pixels of specified pixsize.
    
    npix      - x, y dim of square sky map
    pixsize   - pixel size of sky map (radians)
    wn        - wavenumber of observation (cm-1)
    bxbylist  - list of baselines observed. Each baseline
                stored as a tuple (bx, by) (metres)

    returns:    
    dirty_beam - npix by npix numpy array with dirty beam.
    clean_beam - npix by npix numpy array with clean beam.

    The dirty beam is constructed by calculating the Fourier
    transform of the input list of measured baselines.

    A Gaussian is fitted to the centre of the dirty beam and
    the clean beam constructed from the fitted parameters.
    """

    lamb = 1.0 / (wn * 100.0)

    # maximum baseline in the UV plane corresponds to the cube pixsize
    # at this wavelength
    maxbx = lamb / (2.0 * pixsize)

    dirty_beam = numpy.zeros([npix, npix])

    # get grid indices relative to centre of array [npix/2, npix/2]
    grid = numpy.indices((npix, npix))
    grid -= (npix/2 - 1)

    # build up Fourier transform by coadding cosines corresponding to
    # each baseline in the uv plane
    rpix = npix / 2
    for bxby in bxbylist:

        # Used following to test algorithm. This baseline should
        # produce a cosine wave with 1 cycle covering the x extent
        # of the image.
        #bxby = (lamb / (2.0 * rpix * pixsize), 0.0)

        # length and angle of baseline relative to x axis
        length = numpy.sqrt(bxby[0]**2 + bxby[1]**2) / maxbx
        theta = numpy.arctan2(bxby[1], bxby[0])

        contribution = grid[1] * numpy.cos(theta) + grid[0] * numpy.sin(theta)
        contribution = numpy.cos(contribution * numpy.pi * length)

        dirty_beam += contribution

    # normalise
    dirty_beam /= numpy.max(dirty_beam)

    # debug prints to show centre of beam - should be 1.0 at [npix/2-1,npix/2-1].
    #centre = npix/2-1
    #print dirty_beam[centre-1:centre+2, centre-1:centre+2]

    # fit Gaussian to centre of dirty beam and use this as the 'clean beam'
    shape = numpy.shape(dirty_beam)
    fitter = fitgaussian.FitGaussian()
    dirty_centre = numpy.array(dirty_beam[shape[0]/2-5:shape[0]/2+5,
      shape[1]/2-5:shape[1]/2+5])
    p = fitter.fitgaussian(dirty_centre)

    # construct the clean beam
    cp = (1.0, float(shape[0])/2.0, float(shape[1])/2.0, p[3], p[4], p[5])
    rotgauss = fitter.gaussian(*cp)
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
        print 'Calculating synthesis dirty beam and clean beam...'

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
        rpix = npix / 2

        mx,my = np.meshgrid(np.arange(-rpix, rpix), np.arange(-rpix, rpix))
        self.result['dirty beam'] = collections.OrderedDict()
        self.result['clean beam'] = collections.OrderedDict()

        # calculate synthesis beams using pp
        jobs = {}
        wnmax = np.max(wn)
        for wavenum in wn:
            #
            # Submit pp calls to calculate_synthesis_beams to calculate
            # the results.
            indata = (npix, pixsize, wavenum, uvmapgen['bxby'],)
            jobs[wavenum] = self.job_server.submit(
              calculate_synthesis_beams, indata, (),
              ('numpy','fitgaussian',))

        for wavenum in wn:
            if jobs[wavenum]() is None:
                raise Exception, 'calculate_synthesis_beams has failed'

            dirty_beam,clean_beam = jobs[wavenum]()
            image = co.Image(data=dirty_beam, axes=[axis1, axis2],
              title='Dirty Beam %06.4g cm-1' % wavenum)
            self.result['dirty beam'][wavenum] = image
            image = co.Image(data=clean_beam, axes=[axis1, axis2],
              title='Clean Beam %06.4g cm-1' % wavenum)
            self.result['clean beam'][wavenum] = image

        return self.result

    def __repr__(self):
        return 'SynthesisBeamsGenerator'

