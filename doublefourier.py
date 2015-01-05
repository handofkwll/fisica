from __future__ import absolute_import

import collections
import math
import numpy as np
import pp

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


class DoubleFourier(object):
    """Class to compute interferograms.
    """

    def __init__(self, parameters, previous_results, job_server):
        self.parameters = parameters
        self.previous_results = previous_results
        self.job_server = job_server

        self.result = collections.OrderedDict()      

    def run(self):
        print 'DoubleFourier.run'

        # convert lengths to m. Wavenumbers leave as cm-1.
        fts = self.previous_results['fts']
        fts_wn = fts['fts_wn']
        fts_wn_truncated = fts['fts_wn_truncated']
        opd_max = fts['opd_max'] / 100.0
        fts_nsample = fts['ftsnsample']
        vdrive = fts['vdrive'] / 100.0
        delta_opd = fts['delta_opd']
        delta_time = delta_opd / vdrive

        times = np.arange(int(fts_nsample), dtype=np.float)
        times *= delta_time
        opd_start = -((fts_nsample + 1) / 2) * delta_opd

        beamsgenerator = self.previous_results['beamsgenerator']

        uvmapgenerator = self.previous_results['uvmapgenerator']
        bxby = uvmapgenerator['bxby']

        skygenerator = self.previous_results['skygenerator']
        skymodel = skygenerator['sky model']
        spatial_axis = self.result['spatial axis'] = skygenerator['spatial axis']
        self.result['frequency axis'] = fts_wn_truncated

        # assuming nx is even then transform has 0 freq at origin and [nx/2] is
        # Nyquist frequency. Nyq freq = 0.5 * Nyquist sampling freq. 
        # Assume further that the fft is shifted so that 0 freq is at nx/2
        nx = len(spatial_axis)
        spatial_freq_axis = np.arange(-nx/2, nx/2, dtype=np.float)
        sample_freq = (180.0 * 3600.0 / np.pi) / (spatial_axis[1] - spatial_axis[0])
        spatial_freq_axis *= (sample_freq / nx)
        self.result['spatial frequency axis'] = spatial_freq_axis
        
        self.result['baseline interferograms'] = collections.OrderedDict()
#        for baseline in bxby:
        for baseline in bxby[:1]:
            print 'baseline', baseline
            measurement = np.zeros(np.shape(times))
            opd = np.zeros(np.shape(times))

            # FTS path diff and possibly baseline itself vary with time
            for tindex,t in enumerate(times):

                # calculate the sky that the system is observing at this
                # moment, incorporating various errors in turn. Be explicit
                # about the copy otherwise a reference is made instead.
                sky_now = skymodel.copy()

                # 1. baseline should be perp to centre of field.
                #    If baseline is tilted then origin of sky map shifts.
                #    (I think effect could be corrected by changing
                #    FTS sample position to compensate.?)
                # for now assume 0 error but do full calculation for timing
                # purposes

                # perfect baseline is perpendicular to direction to 'centre'
                # on sky. Add errors in x,y,z - z towards sky 'centre'
                bz = 0.0
                # convert bz to angular error
                blength = math.sqrt(baseline[0]*baseline[0] +
                  baseline[1]*baseline[1])
                bangle = (180.0 * 3600.0 / math.pi) * bz / blength 
                bx_error = bangle * baseline[0] / blength
                by_error = bangle * baseline[1] / blength

                # calculate xpos, ypos in units of pixel - numpy arrays 
                # [row,col]
                nx = len(spatial_axis)
                colpos = float(nx-1) * bx_error / \
                  (spatial_axis[-1] - spatial_axis[0])
                rowpos = float(nx-1) * by_error / \
                  (spatial_axis[-1] - spatial_axis[0])

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

                jobs = {}
                # go through freq planes and shift them
                for iwn,wn in enumerate(fts_wn_truncated):
                    # submit jobs
                    indata = (sky_now[:,:,iwn], shift,)
                    jobs[wn] = self.job_server.submit(fftshift,
                      indata, (), ('numpy',))
                
                for iwn,wn in enumerate(fts_wn_truncated):
                    # collect and store results
                    temp = jobs[wn]()
                    sky_now[:,:,iwn] = temp

                if t == times[0]:
                    # take copy of array to fix a snapshot of it now
                    self.result['sky at time 0'] = sky_now.copy()

                # 2. telescopes should be centred on centre of field
                #    Telescopes collect flux from the 'sky' and pass
                #    it to the FTS beam combiner. In doing this each
                #    telescope multiplies the sky emission by its
                #    amplitude beam response - always real but with
                #    negative areas. Is this correct? Gives right
                #    answer for 'no error' case.

                # multiply sky by amplitude beam 1 * amplitude beam 2
                amp_beam_1 = beamsgenerator['primary amplitude beam'].data
                amp_beam_2 = beamsgenerator['primary amplitude beam'].data

                for iwn,wn in enumerate(fts_wn):
                    # calculate shifted beams here
                    # for now assume no errors and just use beam
                    # calculated earlier
                    pass

                # multiply sky by amplitude beams of 2 antennas
                sky_now *= amp_beam_1 * amp_beam_2 

                if t == times[0]:
                    # take copy of array to fix a snapshot of it now
                    self.result['sky*beams at time 0'] = sky_now.copy()
                    
                # 3. baseline error revisited
                # derive baseline at this time
                #
                # Perhaps baseline should be a continuous function in
                # time, which would allow baselines that intentionally
                # smoothly vary (as in rotating tethered assembly)
                # and errors to be handled by one description.
                #
                # what follows assumes zero error

                baseline_error = 0.0
                baseline_now = baseline + baseline_error

                fft_now = np.zeros(np.shape(sky_now), np.complex)
                spectrum = np.zeros(np.shape(fts_wn), np.complex)
                for iwn,wn in enumerate(fts_wn_truncated):

                    # derive shift needed to place baseline at one of FFT coords
                    # this depends on physical baseline and frequency
                    baseline_now_lambdas = baseline_now * wn * 100.0 

                    # calculate baseline position in units of pixels of FFTed 
                    # sky - numpy arrays [row,col]
                    colpos = float(nx-1) * \
                      float(baseline_now_lambdas[0] - spatial_freq_axis[0]) / \
                      (spatial_freq_axis[-1] - spatial_freq_axis[0])
                    rowpos = float(nx-1) * \
                      float(baseline_now_lambdas[1] - spatial_freq_axis[0]) / \
                      (spatial_freq_axis[-1] - spatial_freq_axis[0])
                    colpos = 0.0
                    rowpos = 0.0

                    # calculate fourier phase shift to move point at [rowpos,colpos] to 
                    # [0,0]
                    shiftx = np.zeros([nx], np.complex)
                    shiftx[:nx/2] = np.arange(nx/2, dtype=np.complex)
                    shiftx[nx/2:] = np.arange(-nx/2, 0, dtype=np.complex)
                    shiftx = np.exp((-2.0j * np.pi * colpos * shiftx) / \
                      float(nx))

                    shifty = np.zeros([nx], np.complex)
                    shifty[:nx/2] = np.arange(nx/2, dtype=np.complex)
                    shifty[nx/2:] = np.arange(-nx/2, 0, dtype=np.complex)
                    shifty = np.exp((-2.0j * np.pi * rowpos * shifty) / float(nx))

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

                    # set amp/phase at this frequency
                    spectrum[wn==fts_wn] = temp[0,0]

                if t == times[0]:
                    self.result['skyfft at time 0'] = fft_now.copy()

                    axis = co.Axis(data=fts_wn, title='wavenumber',
                      units='cm-1')
                    temp = co.Spectrum(data=spectrum, axis=axis,
                      title='Detected spectrum', units='W sr-1 m-2 Hz-1')
                    self.result['skyfft spectrum at time 0'] = temp

                # 3. FTS sampling should be accurate
                #    derive lag due to FTS path difference
                #    0 error for now

                if t == times[0]:
                    # test interferogram
                    # inverse fft of emission spectrum at this point
                    reflected_spectrum = np.zeros([2*(len(spectrum)-1)],
                      np.complex)
                    reflected_spectrum[:len(spectrum)].real = spectrum
                    reflected_spectrum[len(spectrum):].real = spectrum[-2:0:-1]
                    
                    temp = np.fft.ifft(reflected_spectrum)
                    pos = np.fft.fftfreq(len(reflected_spectrum), d=fts_wn[1]-fts_wn[0])

                    # move 0 frequency to centre of array
                    temp = np.fft.fftshift(temp)
                    pos = np.fft.fftshift(pos)

                    axis = co.Axis(data=pos, title='path difference',
                      units='cm')
                    temp = co.Spectrum(data=temp, axis=axis,
                      title='Detected interferogram', units='')

                    self.result['test FTS at time 0'] = temp

                mirror_error = 0.0
                opd[tindex] = opd_start + (vdrive * t + mirror_error)
                opd_ipos = opd[tindex] / delta_opd

                # calculate shift needed to move point at opd to 0
                nfreq = len(fts_wn)
                ndoublefreq = 2 * (nfreq - 1)
                shift = np.zeros([ndoublefreq], dtype=np.complex)
                shift[:nfreq] = np.arange(nfreq, dtype=np.complex)
                shift[nfreq:] = np.arange(nfreq, dtype=np.complex)[-2:0:-1]
                shift = np.exp((-2.0j * np.pi * opd_ipos * shift) / \
                  float(ndoublefreq))

                # reflect spectrum about 0 to give unaliased version
                reflected_spectrum = np.zeros([2*(len(spectrum)-1)], np.complex)
                reflected_spectrum[:len(spectrum)].real = spectrum
                reflected_spectrum[len(spectrum):].real = spectrum[-2:0:-1]

                # apply phase shift and fft
                reflected_spectrum *= shift
                spectrum_fft = np.fft.ifft(reflected_spectrum)
                measurement[tindex] = spectrum_fft[0]

#                if np.array_equal(baseline, bxby[0]):
#                    if 'skyfft check' not in self.result.keys():
#                        self.result['skyfft check'] = {}
#                    temp = co.Spectrum(data=spectrum_fft, 
#                      title='Test interferogram', units='')
#                    self.result['skyfft check'][t] = temp

            axis = co.Axis(data=np.array(opd), title='path difference',
              units='cm')
            temp = co.Spectrum(data=measurement, axis=axis,
              title='Detected interferogram', units='')
            self.result['baseline interferograms'][tuple(baseline)] = \
              temp

        return self.result

    def matlab_transform(self):
        # readers should look at Izumi et al. 2006, Applied Optics, 45, 2576
        # for theoretical background. Names of variables in the code correspond
        # to that work.

        # For now, assume 2 light collectors giving one baseline at a time.
        interferograms = {}

        for baseline in self.baselines:
            interferogram = 0

            # baseline length (cm) and position angle
            bu = baseline[0]
            bv = baseline[1]
            mod_b = np.sqrt(pow(bu,2) + pow(bv,2))
            ang_b = np.arctan2(bv, bu)

            # loop over sky pixels covered by primary beam
            nx = self.sky_s.shape[1]
            ny = self.sky_s.shape[2]
            for j in range(ny):
                for i in range(nx):

                    # inverse fft of emission spectrum at this point
                    temp = np.fft.ifft(self.sky_s[:,j,i])

                    # move 0 frequency to centre of array
                    temp = np.fft.fftshift(temp)

                    # length (radians) and position angle of theta vector
                    mod_theta = np.sqrt(pow(self.sky_x[i],2) + pow(self.sky_y[j],2))
                    ang_theta = np.arctan2(self.sky_y[j], self.sky_x[i])

                    # calculate b.theta (the projection of b on theta)
                    # and the corresponding delay in units of wavelength at 
                    # Nyquist frequency
                    delay = mod_theta * mod_b * np.cos(ang_b - ang_theta) * self.freqs[-1]

                    # sampling is done at twice Nyquist freq so shift transformed
                    # spectrum by 2 * delay samples (approximated to nint)
                    # NOTE factor of 2 discrepency with matlab version! I think
                    # this is because there the variable 'Nyq' is the Nyquist
                    # sampling rate, not the Nyquist frequency. 
                    temp = np.roll(temp, int(round(2.0 * delay)))

                    # want only the real part of the result
                    interferogram += np.real(temp) 

    def __repr__(self):
        return 'DoubleFourier'

