from __future__ import absolute_import

import collections
import numpy as np

#used for debugging
#import matplotlib.pyplot as plt
#import common.commonobjects as co

def calculate_visibility(baseline, wn, sky_plane, spatial_freq_axis):
    """Routine to calculate the visibility for a specified
    baseline for one plane in a sky model.

    Parameters:
    baseline          - (u,v) in metres
    wn                - wavenumber of observation (cm-1)
    sky_plane         - 2d sky image at this wn
    spatial_freq_axis - The spatial frequency axis of the sky plane
                        Fourier transform
    """ 

    # derive shift needed to place baseline at one of FFT
    # coords this depends on physical baseline and frequency
    baseline_lambda = baseline * wn * 100.0

    # calculate baseline position in units of pixels of 
    # FFTed sky - numpy arrays [row,col]
    nx,ny = numpy.shape(sky_plane)
    colpos = float(nx-1) * baseline_lambda[0] / \
      (spatial_freq_axis[-1] - spatial_freq_axis[0])
    rowpos = float(nx-1) * baseline_lambda[1] / \
      (spatial_freq_axis[-1] - spatial_freq_axis[0])

    # calculate fourier phase shift to move point at 
    # [rowpos,colpos] to [0,0]
    shiftx = numpy.zeros([nx], numpy.complex)
    shiftx[:nx/2] = numpy.arange(nx/2, dtype=numpy.complex)
    shiftx[nx/2:] = numpy.arange(-nx/2, 0, dtype=numpy.complex)
    shiftx = numpy.exp((2.0j * numpy.pi * colpos * shiftx) / float(nx))

    shifty = numpy.zeros([nx], numpy.complex)
    shifty[:nx/2] = numpy.arange(nx/2, dtype=numpy.complex)
    shifty[nx/2:] = numpy.arange(-nx/2, 0, dtype=numpy.complex)
    shifty = numpy.exp((2.0j * numpy.pi * rowpos * shifty) / float(nx))

    shift = numpy.ones([nx,nx], numpy.complex)
    for j in range(nx):
        shift[j,:] *= shiftx
    for i in range(nx):
        shift[:,i] *= shifty

    # move centre of sky image to origin
    temp = numpy.fft.fftshift(sky_plane)

    # apply phase shift
    temp *= shift
    # 2d fft
    temp = numpy.fft.fft2(temp)
                  
    # return the visibility
    return temp[0,0]

def fftshift(data, shift):
    """Method to shift a 2d complex data array by applying the
       given phase shift to its Fourier transform.

    Parameters:
    data  - 2d data to be shifted
    shift - 2d array containing phase shift
    """
    # move centre of image to array origin
    temp = numpy.fft.fftshift(data)
    # 2d fft
    temp = numpy.fft.fft2(temp)
    # apply phase shift
    temp *= shift
    # transform and shift back
    temp = numpy.fft.ifft2(temp)
    temp = numpy.fft.fftshift(temp)

    return temp

def fields_match(previous_config, config, exclude_fields=[]):
    """Function to compare 2 configs, ignoring any differences in
    the components listed in 'exclude_fields'.

    Parameters:
    previous_config - namedtuple containing the previous configuration
    config          - namedtuple with the configuration to compare
    exclude_fields  - list of tuple fields to be excluded in the
                      comparison
    """
    if previous_config is None:
        return False

    if not exclude_fields:
        return previous_config==config
    else:
        # convert to dicts and eliminate exclude_fields
        previous_dict = previous_config._asdict()
        current_dict = config._asdict()
        for field in exclude_fields:
            ignore = previous_dict.pop(field, None)
            ignore = current_dict.pop(field, None)

        previous_set = set([k for k in previous_dict.items()])
        current_set = set([k for k in current_dict.items()])
        diff = previous_set.symmetric_difference(current_set)

        debug = False
        if(debug):
            print diff

        return not(bool(diff))


class Observe(object):
    """Class to compute interferograms.
    """

    def __init__(self, parameters, previous_results, job_server):
        """Constructor.

        Parameters:
        parameters       - Dict with parameters from the Excel files.
        previous_results - Current results structure of the simulation run.
        job_server       - ParallelPython job server. 
        """
        self.parameters = parameters
        self.previous_results = previous_results
        self.job_server = job_server

        self.result = collections.OrderedDict()      

    def run(self):
        """Method invoked to calculate the interferograms.
        """
#        print 'Observe.run'

        # access primary beam information
#        primarybeams = self.previous_results['primarybeams']

        # access required FTS information
        fts = self.previous_results['fts']
        fts_wn = fts['fts_wn']
        fts_wn_truncated = fts['fts_wn_truncated']
        # convert delta_opd from cm to m
        delta_opd = fts['delta_opd'] / 100.0
        smec_opd_to_mpd = fts['smec_opd_to_mpd']

        # access to model sky
        skygenerator = self.previous_results['skygenerator']
        sky_model = skygenerator['sky model']
        spatial_axis = self.result['spatial axis'] = skygenerator['spatial axis']
        nx = len(spatial_axis)

        # get list of instrument configuratons       
        timeline = self.previous_results['timeline']
        obs_timeline = timeline['obs_timeline']

        # Calculate the axes of spatial FFTs of the sky model.
        # Assuming nx is even then the transform of a spatial axis has 0 freq 
        # at origin and [nx/2] is Nyquist frequency 
        # [Nyq freq = 0.5 * Nyquist sampling freq].
        # The fft will be shifted so that 0 freq is at nx/2.
        nx = len(spatial_axis)
        spatial_freq_axis = np.arange(-nx/2, nx/2, dtype=np.float)
        # spatial axis in arcsec
        sample_freq = (180.0 * 3600.0 / np.pi) / (spatial_axis[1] - spatial_axis[0])
        spatial_freq_axis *= (sample_freq / nx)
        self.result['spatial frequency axis'] = spatial_freq_axis

        # calculate the measured result for each configuration
        previous_config = None
        observed_times = obs_timeline.keys()
        observed_times.sort()

        debug_plotted = False

        for t in observed_times:
            config = obs_timeline[t]
  
            # Does the spatial configuration for this time match that
            # of the previous one?

# commented out test would calculate a 'new sky' whenever
# the flux collector pointings change. This makes things too
# slow - needs more thought
#            if fields_match(previous_config, config,
#              exclude_fields=['scan_number',
#                              'time',
#                              'baseline_x',
#                              'baseline_y',
#                              'baseline_z',
#                              'baseline_number',
#                              'smec_position',
#                              'smec_nominal_position',
#                              'smec_vel_error',                           
#                              'flag',
#                              'data']):

            if fields_match(previous_config, config,
              exclude_fields=['scan_number',
                              'time',
                              'baseline_x',
                              'baseline_y',
                              'baseline_z',
                              'baseline_number',
                              'pointing1_x',
                              'pointing1_y',
                              'pointing2_x',
                              'pointing2_y',
                              'smec_position',
                              'smec_nominal_position',
                              'smec_vel_error',                           
                              'flag',
                              'data']):

                # yes, then just reuse relevant results from previous
                # config.
                pass

            else:
#                print 'previous', previous_config
#                print 'current', config
                print 'calculating new sky'

                # no, calculate new 'sky'

                # Be explicit about copies otherwise a reference is
                # made which corrupts the original and leads to much
                # confusion.
                sky_now = sky_model.copy()

                # calculate the sky that the system is observing at this
                # moment, incorporating various errors in turn.

                # 1. baseline should be perp to vector towards centre
                #    of field. If baseline is tilted then origin of sky 
                #    map shifts.
                #    I think this error can be folded into that for the
                #    FTS mirror, so leave handling it until then. Is
                #    this correct?
                
                # 2. Flux collectors should be pointing at centre of 
                #    field. They each collect flux from the 'sky' and 
                #    pass it to the FTS beam combiner.
                #    In doing this each telescope multiplies the sky 
                #    emission by its complex amplitude beam response.
                #    So instead of measuring the correlation of the
                #    'sky' with itself, we measure that x Beam1.Beam2*. 
                #    Is this correct? Gives right answer for 'no error'
                #    where Beam1.Beam2* == PSF.

                # multiply sky by amp.beam 1 * conj(amp.beam 2)
#                amp_beam_1 = primarybeams['primary amplitude beam'].data
#                amp_beam_2 = primarybeams['primary amplitude beam'].data
                amp_beam_1 = 1.0
                amp_beam_2 = 1.0

                # shift beams if pointing error
                if config.pointing1_x or config.pointing1_y:
                    # calculate shifted beam1
#                    raise Exception, 'not implemented'
                    current_amp_beam_1 = amp_beam_1
                else:
                    current_amp_beam_1 = amp_beam_1

                if config.pointing2_x or config.pointing2_y:
                    # calculate shifted beam2
#                    raise Exception, 'not implemented'
                    current_amp_beam_2 = amp_beam_2
                else:
                    current_amp_beam_2 = amp_beam_2
 
                sky_now *= current_amp_beam_1 * \
                  np.conj(current_amp_beam_2)            

            # 3. Get result for baseline at this time

# commented out test would calculate a new baseline spectrum
# whenever the collector pointings changed. Too slow - need
# a different approach
#            if fields_match(previous_config, config,
#              exclude_fields=['scan_number',
#                              'time',
#                              'smec_position',
#                              'smec_nominal_position',
#                              'smec_vel_error',                           
#                              'flag',
#                              'data']):

            if fields_match(previous_config, config,
              exclude_fields=['scan_number',
                              'time',
                              'smec_position',
                              'smec_nominal_position',
                              'smec_vel_error',                           
                              'pointing1_x',
                              'pointing1_y',
                              'pointing2_x',
                              'pointing2_y',
                              'flag',
                              'data']):
                pass

            else:
                # calculate new baseline spectrum
                baseline = np.array([config.baseline_x, config.baseline_y])
                print 'calculating new baseline spectrum', baseline
 
                # submit jobs
                jobs = {}
                for iwn,wn in enumerate(fts_wn_truncated):
                    indata = (baseline, wn, sky_now[:,:,iwn],
                      spatial_freq_axis,)
                    jobs[wn] = self.job_server.submit(calculate_visibility,
                      indata, (), ('numpy',))

                # collect results
                baseline_spectrum = np.zeros(np.shape(fts_wn), np.complex)
                for wn in fts_wn_truncated:
                    if jobs[wn]() is None:
                        raise Exception, 'calculate_visibility has failed'

                    # set amp/phase at this frequency
                    baseline_spectrum[wn==fts_wn] = jobs[wn]()

            # 4. Calculate interferogram value from this baseline on the 'sky'
            #    at the current FTS mirror position.

            opd = config.smec_position / smec_opd_to_mpd
            opd_ipos = opd / delta_opd

            # calculate shift needed to move point at opd to 0
            nfreq = len(fts_wn)
            ndoublefreq = 2 * (nfreq - 1)
            shift = np.zeros([ndoublefreq], dtype=np.complex)
            shift[:nfreq] = np.arange(nfreq, dtype=np.complex)
            shift[nfreq:] = np.arange(nfreq, dtype=np.complex)[-2:0:-1]
            shift = np.exp((-2.0j * np.pi * opd_ipos * shift) /\
              float(ndoublefreq))

            # reflect spectrum about 0 to give unaliased version
            reflected_spectrum = np.zeros([2*(len(baseline_spectrum)-1)],
              np.complex)
            reflected_spectrum[:len(baseline_spectrum)] = baseline_spectrum
            reflected_spectrum[len(baseline_spectrum):] = \
              baseline_spectrum[-2:0:-1]

            # apply phase shift and fft
            reflected_spectrum *= shift
            spectrum_fft = np.fft.ifft(reflected_spectrum)
            config = config._replace(data = spectrum_fft[0].real)
            obs_timeline[t] = config

            previous_config = config

        self.result['observed_timeline'] = obs_timeline
        return self.result

    def __repr__(self):
        return '''
Observe:
  timelength length: {timeline_len}
'''.format(
          timeline_len=len(self.result['observed_timeline']))

