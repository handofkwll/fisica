from __future__ import absolute_import

import collections
import math
import numpy as np
import pp

import matplotlib.pyplot as plt

import common.commonobjects as co

def fftshift(data, shift):
    """Method to shift a 2d complex data array by applying the
       given phase shift to its Fourier transform.
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
        self.parameters = parameters
        self.previous_results = previous_results
        self.job_server = job_server

        self.result = collections.OrderedDict()      

    def run(self):
        print 'Observe.run'

        # access primary beam information
        primarybeams = self.previous_results['primarybeams']

        # access required FTS information
        fts = self.previous_results['fts']
        fts_wn = fts['fts_wn']
        fts_wn_truncated = fts['fts_wn_truncated']
        delta_opd = fts['delta_opd']

        # access to model sky
        skygenerator = self.previous_results['skygenerator']
        sky_model = skygenerator['sky model']
        spatial_axis = self.result['spatial axis'] = skygenerator['spatial axis']
        nx = len(spatial_axis)

        # assuming nx is even then transform of spatial axis has 0 freq at origin 
        # and [nx/2] is Nyquist frequency [Nyq freq = 0.5 * Nyquist sampling freq].
        # The fft is shifted so that 0 freq is at nx/2.
        nx = len(spatial_axis)
        spatial_freq_axis = np.arange(-nx/2, nx/2, dtype=np.float)
        sample_freq = (180.0 * 3600.0 / np.pi) / (spatial_axis[1] - spatial_axis[0])
        spatial_freq_axis *= (sample_freq / nx)
        self.result['spatial frequency axis'] = spatial_freq_axis

        # get list of instrument configuratons       
        uvmapgenerator = self.previous_results['uvmapgenerator']
        obs_framework = uvmapgenerator['obs_framework']

        # and calculate the result for each configuration
        previous_config = None
        observed_obs_framework = []

        debug_plotted = False
        for config in obs_framework:
#            if config.fts_start:
#                print 'start FTS scan'
#                print 'baseline', config.baseline_x, config.baseline_y,\
#                  config.baseline_z
  
            if fields_match(previous_config, config,
              exclude_fields=['scan_number',
                              'fts_start',
                              'fts_position',
                              'fts_nominal_position',
                              'baseline_x',
                              'baseline_y',
                              'baseline_z',
                              'baseline_number',
                              'data']):
                pass
            else:
#                print 'previous', previous_config
#                print 'current', config
                # calculate new sky
#                print 'calculating new sky'

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
                
                # 2. Flux-collecting telescopes should be pointing 
                #    at centre of field. They each collect flux from 
                #    the 'sky' and pass it to the FTS beam combiner.
                #    In doing this each telescope multiplies the sky 
                #    emission by its complex amplitude beam response.
                #    So instead of measuring the correlation of the
                #    'sky' with itself, we measure that x Beam1.Beam2*. 
                #    Is this correct? Gives right answer for 'no error'
                #    where Beam1.Beam2* == PSF.

                # multiply sky by amp.beam 1 * conj(amp.beam 2)
                amp_beam_1 = primarybeams['primary amplitude beam'].data
                amp_beam_2 = primarybeams['primary amplitude beam'].data

                # shift beams if pointing error
                if config.point1_x_error or config.point1_y_error:
                    # calculate shifted beam1
#                    raise Exception, 'not implemented'
                    current_amp_beam_1 = amp_beam_1
                else:
                    current_amp_beam_1 = amp_beam_1

                if config.point2_x_error or config.point2_y_error:
                    # calculate shifted beam2
#                    raise Exception, 'not implemented'
                    current_amp_beam_2 = amp_beam_2
                else:
                    current_amp_beam_2 = amp_beam_2
 
                sky_now *= current_amp_beam_1 * \
                  np.conj(current_amp_beam_2)            

            # 3. Get result for baseline at this time
            if fields_match(previous_config, config,
              exclude_fields=['fts_start',
                              'fts_position',
                              'fts_nominal_position',
                              'data']):
                pass
            else:
                # calculate new baseline spectrum
#                print 'calculating new baseline spectrum'
 
                baseline = np.array([config.baseline_x, config.baseline_y])
                fft_now = np.zeros(np.shape(sky_now), np.complex)
                baseline_spectrum = np.zeros(np.shape(fts_wn), np.complex)

                for iwn,wn in enumerate(fts_wn_truncated):

                    # derive shift needed to place baseline at one of FFT
                    # coords this depends on physical baseline and frequency
                    baseline_lambda = baseline * wn * 100.0

                    # calculate baseline position in units of pixels of 
                    # FFTed sky - numpy arrays [row,col]
                    colpos = float(nx-1) * baseline_lambda[0] / \
                      (spatial_freq_axis[-1] - spatial_freq_axis[0])
                    rowpos = float(nx-1) * baseline_lambda[1] / \
                      (spatial_freq_axis[-1] - spatial_freq_axis[0])

                    # calculate fourier phase shift to move point at 
                    # [rowpos,colpos] to [0,0]
                    shiftx = np.zeros([nx], np.complex)
                    shiftx[:nx/2] = np.arange(nx/2, dtype=np.complex)
                    shiftx[nx/2:] = np.arange(-nx/2, 0, dtype=np.complex)
                    shiftx = np.exp((2.0j * np.pi * colpos * shiftx) / \
                      float(nx))

                    shifty = np.zeros([nx], np.complex)
                    shifty[:nx/2] = np.arange(nx/2, dtype=np.complex)
                    shifty[nx/2:] = np.arange(-nx/2, 0, dtype=np.complex)
                    shifty = np.exp((2.0j * np.pi * rowpos * shifty) / \
                      float(nx))

                    shift = np.ones([nx,nx], np.complex)
                    for j in range(nx):
                        shift[j,:] *= shiftx
                    for i in range(nx):
                        shift[:,i] *= shifty

                    # move centre of sky image to origin
                    temp = np.fft.fftshift(sky_now[:,:,iwn])
                    # apply phase shift
                    temp *= shift
                    # 2d fft
                    temp = np.fft.fft2(temp)
                    fft_now[:,:,iwn] = temp

                    if wn==30.0 and not debug_plotted:

                        freq = np.fft.fftfreq(n=np.shape(fft_now)[0],
                          d=(spatial_axis[1]-spatial_axis[0]) / 206265.0)
                        freq = np.arange(np.shape(fft_now)[0]) * (freq[1] - freq[0])

                        ifreq = np.arange(np.shape(fft_now)[0])
                        baseline_i = ifreq[abs(freq-abs(baseline_lambda[0])) < (freq[1]-freq[0])]
                        baseline_j = ifreq[abs(freq-abs(baseline_lambda[1])) < (freq[1]-freq[0])]
#                        print baseline, fft_now[baseline_i[0], baseline_j[0], iwn]

#                        plt.figure()
#                        plt.imshow(fft_now.real[:,:,iwn], interpolation='nearest', origin='lower',
#                          aspect='equal', extent=[freq[0],freq[-1],freq[0],freq[-1]])
#                        plt.colorbar(orientation='vertical')
#                        plt.axis('image')
#                        plt.savefig('debug_fft_real.png')
#                        plt.close()

#                        plt.figure()
#                        plt.imshow(fft_now.imag[:,:,iwn], interpolation='nearest', origin='lower',
#                          aspect='equal', extent=[freq[0],freq[-1],freq[0],freq[-1]])
#                        plt.colorbar(orientation='vertical')
#                        plt.axis('image')
#                        plt.savefig('debug_fft_imag.png')
#                        plt.close()

#                        debug_plotted = True
                  
                    # set amp/phase at this frequency
                    baseline_spectrum[wn==fts_wn] = temp[0,0]

            # 4. Calculate interferogram value from this baseline on the 'sky'
            #    at the current FTS mirror position.

            opd = 2.0 * config.fts_position
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
      
            observed_obs_framework.append(config)

            previous_config = config

        self.result['observed_framework'] = observed_obs_framework
        return self.result

    def __repr__(self):
        return 'Observe'

