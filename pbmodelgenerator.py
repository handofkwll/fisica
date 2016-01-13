"""This module contains classes and methods used to calculate the
primary beams.
"""

from __future__ import absolute_import

import collections
import numpy as np
import os
import re

import common.commonobjects as co

def _getfilelist(filedir, fileroot):
    """Utility routine to return a list of all files in
    directory 'filedir' whose names start with the string in
    'fileroot'.

    Keyword parameters:
    filedir  - The name of the directory containing the files.
    fileroot - Returned files will have names starting with this
               string.

    Returns:
             - A list of all the files in 'filedir' whose
               names begin with the string in 'fileroot'.
    """
    files = os.listdir(filedir)
    regexp = re.compile('%s.*' % fileroot)
    files = [file for file in files if regexp.match(file)]
    return files

def _get_baseline_wavelength_models(filedir, filelist):
    """Utility routine to process a list Maynooth
    beam model files with names
    of form ..._bas<baseline>_wav<wavelength>...
    and return a dictionary of the model data for
    each (baseline, wavelength) present.

    Keyword parameters:
    filedir  - The name of the directory containing the files.
    filelist - List of files names of form 
               <blah>_bas<baseline>_wav<wavelength>

    Returns:
             - A dictionary whose keys are tuples
               (baseline, wavelength) derived by parsing
               the names in 'filelist'. The value at
               each key is the illumination model read in
               from that file. This is itself a dict
               with keys 'limits', 'ex', 'ey' and 'ez'
               holding the extent of the modelled area
               and matrices with the x, y, z components
               of the E field. 
             - rotation resolution [deg]. When the angle of the
               p.b. on the sky changes by more than this value then
               the beam should be 'rotated' to suit.
    """
    result = {}
    baseline_re = re.compile('.*bas([0-9]+).*')
    wave_re = re.compile('.*wav([0-9]+).*')
    for file in filelist:
        try:
            # parse the filename to find the baseline/wavelength
            baseline = baseline_re.match(file).groups()
            baseline = float(baseline[0])
            wavelength = wave_re.match(file).groups()
            wavelength = float(wavelength[0])
        except:
            raise Exception, \
              'error parsing baseline and wavelength from filename: %s'\
              % file 

        assert (baseline, wavelength) not in result.keys()

        # read the model
        xmin, ymin, xmax, ymax, ex, ey, ez = _readpbfile(filedir, file)

        # wavelength in microns, wn in cm-1
        wn = 1.0e4 / wavelength
        result[(baseline, wn)] = {
          'limits': (xmin, ymin, xmax, ymax),
          'ex':ex, 'ey':ey, 'ez':ez}

    # these models are generally not circularly symmetric so  
    # need to be rotated with the baseline. Set the angular shift
    # requiring the beam to recalculated.
    rotation_resolution = 5.0

    return result, rotation_resolution

def _readpbfile(filedir, pbfile):
    """Utility routine to read in a primary beam model from
    a Maynooth format file.

    Keyword parameters:
    filedir  - The name of the directory containing the file.
    pbfile   - The name of the file containing the primary beam
               model.

    Returns:
             - xmin. The minimum of the x extent covered.
             - ymin. The minimum of the y extent covered.
             - xmax. The maximum of the x extent covered.
             - ymax. The maximum of the y extent covered.
             - ex_cube. Cube object containing E_x at each x,y.
               The third cube dim is degenerate. 
             - ey_cube. Cube object containing E_y at each x,y.
             - ez_cube. Cube object containing E_z at each x,y.
    """
    f = open(os.path.join(filedir, pbfile), mode='r')
    line_no = 0
    for line in f:
        #print line, line_no
        if line_no == 5:
            assert '++++' in line
        elif line_no == 9:
            try:
                numbers = line.split()
                xmin = float(numbers[0]) 
                ymin = float(numbers[1]) 
                xmax = float(numbers[2]) 
                ymax = float(numbers[3])
            except:
                f.close()
                raise Exception, 'failed to read xmin, ymin, xmax, ymax: %s' % line
        elif line_no == 10:
            try:
                numbers = line.split()
                nx = int(numbers[0]) 
                ny = int(numbers[1])
                ex = np.zeros([nx, ny], np. complex)
                ravel_ex = np.ravel(ex)
                ey = np.zeros([nx, ny], np. complex)
                ravel_ey = np.ravel(ey)
                ez = np.zeros([nx, ny], np. complex)
                ravel_ez = np.ravel(ez)
                ravel_index = 0
            except:
                f.close()
                raise Exception, 'failed to read nx, ny: %s' % line
        elif line_no > 10:
            try:
                numbers = line.split()
                ravel_ex[ravel_index] = complex(float(numbers[0]),
                  float(numbers[1])) 
                ravel_ey[ravel_index] = complex(float(numbers[2]),
                  float(numbers[3])) 
                ravel_ez[ravel_index] = complex(float(numbers[4]),
                  float(numbers[5])) 
                ravel_index += 1
            except:
                f.close()
                raise Exception, 'failed to decode data line: %s' % line

        line_no += 1

    # store the model results in Cube objects
    axis = np.arange(nx, dtype=np.float)
    assert nx%2 == 1
    axis = axis - ((nx-1) / 2)
    axis *= (xmax - xmin) / float(nx-1)
    axis1 = co.Axis(data=axis, title='x', units='m') 
    axis2 = co.Axis(data=axis, title='y', units='m') 

    ex_cube = co.Cube(data=ex, axes=[axis1, axis2], title='$E_x$')
    ey_cube = co.Cube(data=ey, axes=[axis1, axis2], title='$E_y$')
    ez_cube = co.Cube(data=ez, axes=[axis1, axis2], title='$E_z$')
  
    f.close()

    return xmin, ymin, xmax, ymax, ex_cube, ey_cube, ez_cube

def _get_perfect_baseline_wavelength_model(m1_diameter):
    """Utility routine to construct a Maynooth type
    result that is just even illumination across the
    collector primary.

    Keyword parameters:
    m1_diameter - Diameter of collector primary in metres.

    Returns:
             - A dictionary with one key (50.0, 100.0)
               whose value is a model with even 
               illumination across the collector 
               primary.
             - rotation resolution [deg]. When the angle of the
               p.b. on the sky changes by more than this value then
               the beam should be 'rotated' to suit.
    """
    result = {}

    # extent similar to Maynooth models
    xmin = -0.6 * m1_diameter
    ymin = -0.6 * m1_diameter
    xmax = 0.6 * m1_diameter
    ymax = 0.6 * m1_diameter
    nx = 101

    # the 'perfect' illumination model is a circle
    # of ones over the mirror area, zeros outside
    exy = np.zeros([nx,nx], np.complex)
    pos = np.indices([nx,nx], np.float)
    pos -= (nx-1)/2
    pos *= (xmax - xmin) / (nx-1)
    exy[(pos[0]**2 + pos[1]**2) < (m1_diameter / 2.0)**2] = 1.0
    ez = np.zeros([nx,nx], np.complex)

    # store the results in Cube objects
    axis = np.arange(nx, dtype=np.float)
    assert nx%2 == 1
    axis = axis - ((nx-1) / 2)
    axis *= (xmax - xmin) / float(nx-1)
    axis1 = co.Axis(data=axis, title='x', units='m') 
    axis2 = co.Axis(data=axis, title='y', units='m') 

    ex_cube = co.Cube(data=exy, axes=[axis1, axis2], title='$E_x$')
    ey_cube = co.Cube(data=exy, axes=[axis1, axis2], title='$E_y$')
    ez_cube = co.Cube(data=ez, axes=[axis1, axis2], title='$E_z$')

    models = {}
    models[(50.0, 100.0)] = {
      'limits': (xmin, ymin, xmax, ymax),
      'ex':ex_cube, 'ey':ey_cube, 'ez':ez_cube}

    # these models are circularly symmetric so do not need to be
    # rotated with the baseline. Set the angular shift requiring
    # the beam be rotated to 360 degrees (i.e. never triggered).
    rotation_resolution = 360.0

    return models, rotation_resolution

def calculate_primary_beam_from_pbmodel(npix, pixsize, m1_diameter, wn,
  pbmodel, xmin, ymin, xmax, ymax,):
    """Routine to calculate the primary beams on a 
    sky map with npix x npix pixels of specified pixsize.
    
    Keyword parameters:
    npix        - x, y dim of square sky map
    pixsize     - pixel size of sky map (radians)
    m1_diameter - Diameter of the flux collector primary
                  mirrors (metres).
    wn          - wavenumber of observation (cm-1).
    pbmodel     - complex E field just in front of primary.
    xmin        - minimum x coord of pbmodel.
    ymin        - minimum y coord of pbmodel.
    xmax        - maximum x coord of pbmodel.
    ymax        - maximum y coord of pbmodel.

    Returns:    
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

    ny,nx = numpy.shape(pbmodel)
    assert nx%2 == 1
    ix0 = (nx-1) / 2
    iy0 = (ny-1) / 2
    dx = (xmax - xmin) / float(nx-1)
    dy = (ymax - ymin) / float(ny-1)

    # Calculate the Fourier transform.
    # .. pbmodel stored [x,y], y varies fastest
    # .. only use every 3rd index to speed things up by a factor 9
    # .. but take care to use indices symmetric about middle
    step = 3
    start = ((nx-1) / 2) % step
    for ix in numpy.arange(start=start, stop=nx, step=step):
        x = (ix - ix0) * dx
        for iy in numpy.arange(start=start, stop=ny, step=step):
            y = (iy - iy0) * dy

            # Used following to test algorithm. This baseline should
            # produce a cosine wave with 1 cycle covering the x extent
            # of the image.
            # x = lamb / (2.0 * rpix * pixsize)
            # y = 0.0

            # length and angle of baseline relative to x axis
            length = numpy.sqrt(x**2 + y**2) / maxbx
            theta = numpy.arctan2(y, x)

            contribution = grid[1] * numpy.cos(theta) + grid[0] * numpy.sin(theta)
            contribution = numpy.exp(1j * contribution * numpy.pi * length)
            contribution *= pbmodel[iy,ix]

            primary_amplitude_beam += contribution

    # normalise
    beam_max = numpy.max(numpy.abs(primary_amplitude_beam))
    if beam_max > 0.0:
        primary_amplitude_beam /= beam_max

    return primary_amplitude_beam


class PrimaryBeamsGenerator(object):
    """Class to generate the primary beam(s) of the simulated observation.

    Contains methods:
    __init__
    _calculate_beam_
    run
    __repr__
    """

    def __init__(self, previous_results, beam_model_dir, job_server):
        """Constructor.

        Keyword parameters:
        previous_results - Results structure of the simulation run.
        beam_model_dir   - Directory to be searched for beam models.
        job_server       - ParallelPython job server.
        """
        self.previous_results = previous_results
        self.job_server = job_server

        self.beam_model_dir = beam_model_dir

        self.result = collections.OrderedDict()
        self.result['beam_model_dir'] = beam_model_dir

    def _calculate_beam(self, beam_model_type, beam_model_dir, 
      beam_model_pol, wn, npix, pixsize, m1_diameter):
        """Private method invoked to calculate the beams.

        Keyword parameters:
        beam_model_type - 'perfect' will give the beams for a perfectly
                          illuminated primary, any other string is the rootname
                          of a set of beam models in beam_model_dir that are
                          to be used.
        beam_model_dir  - directory to be searched for beam models
        beam_model_pol  - polarization of model to be used: 'X', 'Y' or 'Z'
        wn              - frequencies for which beams to be calculated [cm-1]
        npix            - beam maps have dimensions [npix,npix]
        pixsize         - pixel size of beam maps [rad]
        m1_diameter     - diameter of primary mirror [m]

        Returns:
        intensity_beams - intensity beam (Amp.Amp*) cube [nwn, npix, npix]
        amplitude_beams - complex amplitude beam cube [nwn, npix, npix]
        models          - dictionary containing beam models used
        rotation_res    - [deg]. When the angle of the
                          p.b. on the sky changes by more than this value then
                          the beam should be 'rotated' to suit.
        """

        # spatial axes same for all wavelengths/model types
        rpix = npix / 2
        axis = np.arange(-rpix, rpix, dtype=np.float)
        axis *= pixsize
        axis = np.rad2deg(axis) * 3600.0
        axis1 = co.Axis(data=wn, title='Frequency', units='cm-1')
        axis2 = co.Axis(data=axis, title='Dec offset', units='arcsec') 
        axis3 = co.Axis(data=-axis, title='RA offset', units='arcsec') 

        # containers for results
        intensity_beams = collections.OrderedDict()
        amplitude_beams = collections.OrderedDict()

        if beam_model_type.lower().strip() == 'perfect':
            # in this case the model is not a function of baseline or
            # wavenumber so just calculate a canonical result
            models, rotation_res = _get_perfect_baseline_wavelength_model(
              m1_diameter)

        else:
            # list all files with specified root
            filelist = _getfilelist(beam_model_dir, beam_model_type)

            # get the grid of models available for baselines/wavelengths
            models, rotation_res = _get_baseline_wavelength_models(
              beam_model_dir, filelist)

        primary_illumination = models

        # get list of baseline lengths modelled
        model_grid = models.keys()
        baselines = set()
        for k in model_grid:
            baselines.update([k[0]])
        baselines = list(baselines)
        baselines.sort()

        # get wavenumbers for which models available on this baseline
        wnlist = set()
        for k in model_grid:
            wnlist.update([k[1]])
        wnlist = list(wnlist)
        wnlist.sort()
  
        # calculate beams for each wn on each baseline modelled
        jobs = {}

        # iterate through baseline lengths modelled
        for baseline in baselines:
            print 'baseline', baseline

            # now compute beams for each wn required
            for wavenum in wn:
                #print 'wavenum', wavenum

                # use the illumination model that is closest 
                # in wavenumber
                model_wn = np.array(wnlist)
                model_wn = model_wn[np.argmin(model_wn - wavenum)]

                chosen_model = models[(baseline, model_wn)]
                xmin, ymin, xmax, ymax = chosen_model['limits']
                if beam_model_pol.upper().strip() == 'X':
                    efield = chosen_model['ex'].data
                elif beam_model_pol.upper().strip() == 'Y':
                    efield = chosen_model['ey'].data
                elif beam_model_pol.upper().strip() == 'Z':
                    efield = chosen_model['ez'].data
                else:
                    raise Exception, 'bad model pol: %s' % beam_model_pol

                # submit jobs
                indata = (npix, pixsize, m1_diameter, wavenum, efield, 
                  xmin, ymin, xmax, ymax,)
                jobs[wavenum] = self.job_server.submit(
                  calculate_primary_beam_from_pbmodel,
                  indata, (), ('numpy', 'math',))

            # collect results
            intensity_beam = np.zeros([len(wn),npix,npix], np.float)
            amplitude_beam = np.zeros([len(wn),npix,npix], np.complex)
#            print 'beams set to 1'
#            intensity_beam = np.ones([len(wn),npix,npix], np.float)
#            amplitude_beam = np.ones([len(wn),npix,npix], np.complex)
            for iwn,wavenum in enumerate(wn):
                if jobs[wavenum]() is None:
                    raise Exception, \
                      'calculate_primary_beam_from_pbmodel has failed'

                amplitude_beam[iwn,:,:] = temp = jobs[wavenum]()
                intensity_beam[iwn,:,:] = (temp * np.conjugate(temp)).real

            intensity_beams[baseline] = co.Cube(data=intensity_beam,
              axes=[axis1, axis2, axis3], title='Intensity Beam')
            amplitude_beams[baseline] = co.Cube(data=amplitude_beam,
              axes=[axis1, axis2, axis3], title='Amplitude Beam')

        return intensity_beams, amplitude_beams, models, rotation_res

    def run(self):
        """Method that does the work.
        """
        print 'Calculating primary beams...'

        # collector parameters
        telescope = self.previous_results['loadparameters']['substages']\
          ['Telescope']
        m1_diameter = telescope['Primary mirror diameter']
        self.result['c1_beam_model_type'] = c1_beam_model_type = \
          telescope['Collector 1 beam model type']
        self.result['c2_beam_model_type'] = c2_beam_model_type = \
          telescope['Collector 2 beam model type']
        self.result['beam_model_pol'] = beam_model_pol = \
          telescope['Beam model polarization']

        # gather relevant instrument configuration 
        cubeparams = self.previous_results['skymodel']
        self.result['frequency axis'] = wn = cubeparams['frequency axis']
        self.result['pixsize [rad]'] = pixsize = cubeparams['pixsize [rad]']
        self.result['npix'] = npix = len(cubeparams['spatial axis [arcsec]'])

        # get result for collector 1
        self.result['collector 1 intensity beam'],\
          self.result['collector 1 amplitude beam'],\
          self.result['collector 1 primary illumination'],\
          self.result['collector 1 rotation resolution'] = \
          self._calculate_beam(c1_beam_model_type, self.beam_model_dir,
            beam_model_pol, wn, npix, pixsize, m1_diameter)

        # get result for collector 2
        if c2_beam_model_type == c1_beam_model_type:
            self.result['collector 2 intensity beam'] = \
              self.result['collector 1 intensity beam'] 
            self.result['collector 2 amplitude beam'] = \
              self.result['collector 1 amplitude beam']
            self.result['collector 2 primary illumination'] = \
              self.result['collector 1 primary illumination']
            self.result['collector 2 rotation resolution'] = \
              self.result['collector 1 rotation resolution']
        else:
            self.result['collector 2 intensity beam'],\
              self.result['collector 2 amplitude beam'],\
              self.result['collector 2 primary illumination'],\
              self.result['collector 2 rotation resolution'] = \
              self._calculate_beam(c2_beam_model_type, self.beam_model_dir,
                wn, npix, pixsize, m1_diameter)

        return self.result

    def __repr__(self):
        blurb = '''
PrimaryBeamsGenerator:
  Collector 1:'''

        if self.result['c1_beam_model_type'].lower().strip() == 'perfect':
            blurb += '''
    Evenly illuminated primary'''

        else:
            blurb += '''
    Models read from directory - '{dir}'
    from files with root - '{root}'
    with polarization - '{pol}'
    Baseline Wavenumber'''.format(
              dir=self.result['beam_model_dir'],
              root=self.result['c1_beam_model_type'],
              pol=self.result['beam_model_pol'])

            keys = self.result['collector 1 primary illumination'].keys()
            if keys:
                for k in self.result['collector 1 primary illumination'].keys():
                    blurb += '''
    {baseline}   {wn}'''.format(baseline=k[0],
                              wn=k[1])
            else:
                blurb += '''
    no data read'''

        blurb += '''
    Rotation resolution - {rotres}'''.format(
          rotres=self.result['collector 1 rotation resolution'])

        blurb += '''

  Collector 2:'''

        if self.result['c2_beam_model_type'].lower().strip() == 'perfect':
            blurb += '''
    Evenly illuminated primary'''

        else:
            blurb += '''
    Models read from directory - '{dir}'
    from files with root - '{root}'
    with polarization - '{pol}'
    Baseline Wavenumber'''.format(
              dir=self.result['beam_model_dir'],
              root=self.result['c2_beam_model_type'],
              pol=self.result['beam_model_pol'])

            keys = self.result['collector 2 primary illumination'].keys()
            if keys:
                for k in self.result['collector 2 primary illumination'].keys():
                    blurb += '''
    {baseline}   {wn}'''.format(baseline=k[0],
                              wn=k[1])
            else:
                blurb += '''
    no data read'''

        blurb += '''
    Rotation resolution - {rotres}'''.format(
          rotres=self.result['collector 2 rotation resolution'])

        return blurb
