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

    filedir  - The name of the directory containing the files.
    fileroot - Returned files will have names starting with this
               string.

    returns:
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

    filedir  - The name of the directory containing the files.
    filelist - List of files names of form 
               <blah>_bas<baseline>_wav<wavelength>

    returns:
             - A dictionary whose keys are tuples
               (baseline, wavelength) derived by parsing
               the names in 'filelist'. The value at
               each key is the illumination model read in
               from that file. This is itself a dict
               with keys 'limits', 'ex', 'ey' and 'ez'
               holding the extent of the modelled area
               and matrices with the x, y, z components
               of the E field. 
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
    return result

def _readpbfile(filedir, pbfile):
    """Utility routine to read in a primary beam model from
    a Maynooth format file.

    filedir  - The name of the directory containing the file.
    pbfile   - The name of the file containing the primary beam
               model.

    returns:
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
#    assert ravel_index == nx * ny

    return xmin, ymin, xmax, ymax, ex_cube, ey_cube, ez_cube

def calculate_primary_beam(npix, pixsize, m1_diameter, wn, nuv):
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
        contribution = numpy.exp(1j * contribution * numpy.pi * length)

        primary_amplitude_beam += contribution

        # baselines are symmetric in x so coadd the x 'mirror' of
        # 'contribution' as well.

        primary_amplitude_beam += contribution[::-1]

    # normalise 
    primary_amplitude_beam /= numpy.max(primary_amplitude_beam.real)

    return primary_amplitude_beam

def calculate_primary_beam_from_pbmodel(npix, pixsize, m1_diameter, wn,
  pbmodel, xmin, ymin, xmax, ymax,):
    """Routine to calculate the primary beams on a 
    sky map with npix x npix pixels of specified pixsize.
    
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

    nx,ny = numpy.shape(pbmodel)
    assert nx%2 == 1
    ix0 = (nx-1) / 2
    iy0 = (ny-1) / 2
    dx = (xmax - xmin) / float(nx-1)
    dy = (ymax - ymin) / float(ny-1)

    # Calculate the Fourier transform.
    # .. pbmodel stored [x,y], y varies fastest
    for ix in numpy.arange(nx, step=3):
        x = (ix - ix0) * dx
        for iy in numpy.arange(ny, step=3):
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
            contribution *= pbmodel[ix,iy]

            primary_amplitude_beam += contribution

    return primary_amplitude_beam


class PrimaryBeamsGenerator(object):
    """Class to generate the primary beam(s) of the simulated observation.
    """

    def __init__(self, previous_results, beam_model_dir, job_server):
        """Constructor.
        previous_results - Current results structure of the simulation run.
        job_server       - ParallelPython job server.
        """
        self.previous_results = previous_results
        self.beam_model_dir = beam_model_dir
        self.result = collections.OrderedDict()
        self.job_server = job_server

        self.nuv = 15

    def run(self):
        """Method that does the work.
        """
        print 'Calculating primary beams...'

        # collector parameters
        telescope = self.previous_results['telescope']
        self.result['beam_model_dir'] = beam_model_dir = self.beam_model_dir
        #'/Users/jfl/Dropbox/FP7-FISICA/Topical Telecon Software/Software_Data/Smooth_Walled_Horn/Band 4 GRASP'
        self.result['beam_model_type'] = beam_model_type = telescope['beam_model_type']

        # list all files with specified root
        filelist = _getfilelist(beam_model_dir, beam_model_type)

        # get the grid of models available for baselines/wavelengths
        models = _get_baseline_wavelength_models(beam_model_dir, filelist)
        self.result['primary illumination'] = models

        # gather relevant instrument configuration 
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

        # containers for results
        self.result['primary intensity beam'] = collections.OrderedDict()
        self.result['primary amplitude beam'] = collections.OrderedDict()

        # calculate beams for each wn on each baseline modelled
        jobs = {}

        # iterate through baselines modelled
        model_grid = models.keys()
        baselines = set()
        for k in model_grid:
            baselines.update([k[0]])
        baselines = list(baselines)
        baselines.sort()
  
        for baseline in baselines:
            print 'baseline', baseline

            # get wavenumbers for which models available on this baseline
            wnlist = set()
            for k in model_grid:
                wnlist.update([k[1]])
            wnlist = list(wnlist)
            wnlist.sort()

            # now compute beams for each wn required
            for wavenum in wn[:4]:
#            for wavenum in wn:
                #print 'wavenum', wavenum

                # use the illumination model that is closest 
                # in wavenumber
                model_wn = np.array(wnlist)
                model_wn = model_wn[np.argmin(model_wn - wavenum)]

                chosen_model = models[(baseline, model_wn)]
                xmin, ymin, xmax, ymax = chosen_model['limits']
                ey = chosen_model['ey'].data

                # submit jobs
                indata = (npix, pixsize, m1_diameter, wavenum, ey, 
                  xmin, ymin, xmax, ymax,)
                jobs[wavenum] = self.job_server.submit(
                  calculate_primary_beam_from_pbmodel,
                  indata, (), ('numpy', 'math', 'zernike',))

            # collect results

            primary_intensity_beam = np.zeros([npix,npix,len(wn)],
              np.float)
            primary_amplitude_beam = np.zeros([npix,npix,len(wn)],
              np.complex)
            for iwn,wavenum in enumerate(wn[:4]):
#            for iwn,wavenum in enumerate(wn):
                if jobs[wavenum]() is None:
                    raise Exception, 'calculate_primary_beams has failed'

                primary_amplitude_beam[:,:,iwn] = temp = jobs[wavenum]()
                primary_intensity_beam[:,:,iwn] = (temp * np.conjugate(temp)).real

            cube = co.Cube(data=primary_intensity_beam, axes=[axis1, axis2, axis3],
              title='Intensity Primary Beam')
            self.result['primary intensity beam'][baseline] = cube
            cube = co.Cube(data=primary_amplitude_beam, axes=[axis1, axis2, axis3],
              title='Amplitude Primary Beam')
            self.result['primary amplitude beam'][baseline] = cube

        return self.result

    def __repr__(self):
        blurb = '''
PrimaryBeamsGenerator:
  Models read from directory - '{dir}'
  from files with root - '{root}'
  Baseline Wavenumber'''.format(
          dir=self.result['beam_model_dir'],
          root=self.result['beam_model_type'])

        keys = self.result['primary illumination'].keys()
        if keys:
            for k in self.result['primary illumination'].keys():
                blurb += '''
  {baseline}   {wn}'''.format(baseline=k[0],
                              wn=k[1])
        else:
            blurb += '''
  no data read'''

        return blurb
