from __future__ import absolute_import

import collections
import numpy as np

import common.commonobjects as co

def calculate_primary_beam(npix, pixsize, m1_diameter, wn,):
    # maximum baseline in the UV plane corresponds to the cube pixsize
    # at this wavelength
    lamb = 1.0 / (wn * 100.0)
    bmax = lamb / (2.0 * pixsize)

    # Oversample UV plane by factor 'mult' to reduce aliasing.
    mult = 5

    # Fill uv plane appropriate to primary size, assume
    # uniform illumination
    rpix = npix / 2
    mirror = numpy.zeros([mult*npix, mult*npix])
    phase = numpy.zeros([mult*npix, mult*npix])
    radius = int((0.5 * m1_diameter / bmax) * mult * rpix) + 1
    # numpy array -ve indexing convenient when mirror centred at array origin.
    # Limiting range of of mx and my for efficiency.
    for mx in range(-radius-5, radius+5):
        for my in range(-radius-5,radius+5):
            if mx*mx + my*my < radius * radius:
                rho = math.sqrt(mx*mx + my*my) / radius
                phi = math.atan2(float(my), float(mx))
                mirror[mx,my] = 1.0
#                phase[mx,my] = zernike.zernike(0, 2, rho, phi, norm=False)

    # do a 2d ftt to get primary beam pattern
    temp = numpy.fft.fft2(mirror)
#    temp = phase

    # truncate primary beam pattern and shift to array centre
    primary_amplitude_beam = numpy.zeros([npix, npix], numpy.complex)
    for mx in range(-rpix, rpix):
        for my in range(-rpix,rpix):
            primary_amplitude_beam[rpix+mx,rpix+my] = temp[mx,my]

    # normalise 
    primary_amplitude_beam /= numpy.max(primary_amplitude_beam)

    return primary_amplitude_beam

def zernike1(i, j, m, n):
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

    def run(self):
        print 'Calculating primary beams...'

        # gather configuration 
        cubeparams = self.previous_results['cubeparameters']
        self.result['wn'] = wn = cubeparams['wn']
        self.result['pixsize [rad]'] = pixsize = cubeparams['pixsize [rad]']
        self.result['npix'] = npix = cubeparams['npix']

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
            indata = (npix, pixsize, m1_diameter, wavenum,)
            jobs[wavenum] = self.job_server.submit(calculate_primary_beam,
              indata, (), ('numpy', 'math', 'zernike',))

        # collect results
        self.result['primary beam'] = collections.OrderedDict()
        self.result['primary amplitude beam'] = collections.OrderedDict()

        primary_amplitude_beam = np.zeros([npix,npix,len(wn)],
          np.complex)
        for iwn,wavenum in enumerate(wn):
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
        return 'PrimaryBeamsGenerator'

