from __future__ import absolute_import

import collections
import numpy as np

import common.commonobjects as co

def calculate_primary_beam(npix, pixsize, m1_diameter, wn, nuv,):
    """Routine to calculate the primary beams on a 
    sky map with npix x npix pixels of specified pixsize.
    
    npix        - x, y dim of square sky map
    pixsize     - pixel size of sky map (radians)
    m1_diameter - Diameter of the flux collector primary
                  mirrors (metres).
    wn          - wavenumber of observation (cm-1)
    nuv         - the number of pixels per mirror radius used to 'sample'
                  the uv plane.

    returns:    
    primary_beam - npix by npix numpy complex array with amplitude
                   primary beam.

    The primary beam is constructed by calculating the Fourier
    transform of the primary mirror amplitude/phase profile.
    """

    lamb = 1.0 / (wn * 100.0)

    # maximum baseline in the UV plane corresponds to the cube pixsize
    # at this wavelength
    maxbx = lamb / (2.0 * pixsize)

    primary_amplitude_beam = numpy.zeros([npix, npix], dtype=numpy.complex)

    # get grid indices relative to centre of array [npix/2-1, npix/2-1]
    grid = numpy.indices((npix, npix))
    grid -= (npix/2 - 1)

    # build up Fourier transform by coadding cosines corresponding to
    # each baseline in the uv plane
    rpix = npix / 2

    # uv radius in metres
    radius = 0.5 * m1_diameter

    # construct a list of 'baselines' covered by the primary mirror.
    # Only baselines in the first quadrant - others are included by
    # symmetry.
    bxbylist = []
    for mx in range(nuv):
        bx = radius * float(mx) / float(nuv - 1)
        for my in range(nuv):
            by = radius * float(my) / float(nuv - 1)
            if bx**2 + by**2 < radius**2:
                bxbylist.append((bx, by))

    # Calculate the Fourier transform.
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

        primary_amplitude_beam += contribution

        # baselines are symmetric in x so coadd the x 'mirror' of
        # 'contribution' as well.

        primary_amplitude_beam += contribution[::-1]

    # normalise 
    primary_amplitude_beam /= numpy.max(primary_amplitude_beam.real)

    return primary_amplitude_beam

def zernike1(i, j, m, n):
    """Routine to calculate a zernike polunomial. NOT USED YET.
    """
    ii = np.ravel(i)
    jj = np.ravel(j)
    for k in range(len(ii)):
        print 'zernike', ii[k], jj[k], m, n


class PrimaryBeamsGenerator(object):
    """Class to generate the primary beam(s) of the simulated observation.
    """

    def __init__(self, previous_results, job_server):
        self.previous_results = previous_results
        self.result = collections.OrderedDict()
        self.job_server = job_server

        self.nuv = 15

    def run(self):
        print 'Calculating primary beams...'

        # gather configuration 
        cubeparams = self.previous_results['skymodel']
        self.result['frequency axis'] = wn = cubeparams['frequency axis']
        self.result['pixsize [rad]'] = pixsize = cubeparams['pixsize [rad]']
        self.result['npix'] = npix = len(cubeparams['spatial axis [arcsec]'])

        telescope = self.previous_results['loadparameters']['substages']\
          ['Telescope']
        m1_diameter = telescope['Primary mirror diameter']

        rpix = npix / 2
        # spatial axes same for all wavelengths
        axis = np.arange(-rpix, rpix, dtype=np.float)
        axis *= pixsize
        axis = np.rad2deg(axis) * 3600.0
        axis1 = co.Axis(data=-axis, title='RA offset', units='arcsec') 
        axis2 = co.Axis(data=axis, title='Dec offset', units='arcsec') 
        axis3 = co.Axis(data=wn, title='Frequency', units='cm-1')

        # calculate beams for each point on 'wn'
        jobs = {}

        for wavenum in wn:
            # submit jobs
            indata = (npix, pixsize, m1_diameter, wavenum, self.nuv,)
            jobs[wavenum] = self.job_server.submit(calculate_primary_beam,
              indata, (), ('numpy', 'math', 'zernike',))

        # collect results
        self.result['primary beam'] = collections.OrderedDict()
        self.result['primary amplitude beam'] = collections.OrderedDict()

        primary_amplitude_beam = np.zeros([npix,npix,len(wn)],
          np.complex)
        for iwn,wavenum in enumerate(wn):
            if jobs[wavenum]() is None:
                raise Exception, 'calculate_primary_beams has failed'

            primary_amplitude_beam[:,:,iwn] = temp = jobs[wavenum]()
            primary_intensity_beam = (temp * np.conjugate(temp)).real
            image = co.Image(data=primary_intensity_beam, axes=[axis1, axis2],
              title='Primary Beam %06.4g cm-1' % wavenum)
            self.result['primary beam'][wavenum] = image

        cube = co.Cube(data=primary_amplitude_beam, axes=[axis1, axis2, axis3],
          title='Amplitude Primary Beam')
        self.result['primary amplitude beam'] = cube

        return self.result

    def __repr__(self):

        return '''
PrimaryBeamsGenerator:
  nuv : {nuv}
'''.format(nuv=self.nuv)

