from __future__ import absolute_import

import numpy as np
import os.path
import pyfits

import matplotlib.pyplot as plt


class ZeroErrors(object):
    """Class to generate SMEC positions with no errors.
    """

    def __init__(self, smec_start_pos, smec_velocity, nsamples,
      scan_duration, prescan_duration):
        """Constructor.
        """

        self.smec_start_pos = smec_start_pos
#        self.smec_velocity = smec_velocity
        self.smec_velocity = 0.0
        self.nsamples = nsamples
        self.scan_duration = scan_duration
        self.prescan_duration = prescan_duration

        # construct sampling times
        sample_time  = scan_duration / float(nsamples - 1)
        total_time = prescan_duration + scan_duration
        self.total_nsamples = int(total_time / sample_time) + 1
        self.prescan_nsamples = self.total_nsamples - nsamples
        self.scan_relative_time = np.arange(self.total_nsamples, dtype=np.float)
        self.scan_relative_time *= sample_time

    def _calculate_smec_position(self, scan, scan_start_time, smec_start_pos, 
      direction, relative_time, demand_smec_position, actual_smec_position,
      smec_integrate, vel_error=None):

        if vel_error is None:
            vel_error = np.zeros(np.shape(self.scan_relative_time))

        relative_time[scan*self.total_nsamples : (scan+1)*self.total_nsamples] = \
          scan_start_time + self.scan_relative_time
        # prescan position for now just smec_start_pos
        offset = scan*self.total_nsamples
        demand_smec_position[offset:offset+self.prescan_nsamples] = smec_start_pos
        actual_smec_position[offset:offset+self.prescan_nsamples] = smec_start_pos

        demand_smec_position[offset+self.prescan_nsamples:offset+self.total_nsamples] = \
          smec_start_pos + (self.scan_relative_time[self.prescan_nsamples:] -
          self.scan_relative_time[self.prescan_nsamples]) * \
          self.smec_velocity * direction
#        actual_smec_position[offset+self.prescan_nsamples:offset+self.total_nsamples] = \
#          smec_start_pos + (self.scan_relative_time[self.prescan_nsamples:] -
#          self.scan_relative_time[self.prescan_nsamples]) * \
#          (self.smec_velocity + vel_error[self.prescan_nsamples:]) * direction

        actual_smec_position[offset:offset+self.prescan_nsamples] = 0.0
        actual_smec_position[offset+self.prescan_nsamples:offset+self.total_nsamples] = \
          vel_error[self.prescan_nsamples:]

        smec_integrate[offset+self.prescan_nsamples:offset+self.total_nsamples] = \
          True

    def run(self, nscans):
#        print 'ZeroSmecErrors.run'
        if nscans%2 != 0:
            raise Exception, 'nscans not multiple of 2'

        times = np.zeros([nscans * self.total_nsamples])
        demand_smec_position = np.zeros([nscans * self.total_nsamples])
        actual_smec_position = np.zeros([nscans * self.total_nsamples])
        smec_integrate = np.zeros([nscans * self.total_nsamples], np.bool)

        # first scan will be forward
        direction = 1.0
        scan_start_time = 0.0
        smec_start_pos = self.smec_start_pos

        for scan in range(nscans):

            self._calculate_smec_position(scan, scan_start_time, smec_start_pos, 
              direction, times, demand_smec_position, actual_smec_position,
              smec_integrate)

            # swap direction for next scan
            direction = -direction
            scan_start_time = times[(scan+1)*self.total_nsamples-1]
            smec_start_pos = demand_smec_position[(scan+1)*self.total_nsamples-1]

        return times, demand_smec_position, actual_smec_position, smec_integrate

    def __repr__(self):
        return 'smecposition.ZeroErrors'


class HerschelErrors(ZeroErrors):
    """Class to create SMEC error sequence derived from
    Herschel FTS.
    """

    def __init__(self, smec_start_pos, smec_velocity, nsamples,
      scan_duration, prescan_duration, seed=57,
      jitterspectrum='VelocityJitter1342270191.fits'):

        ZeroErrors.__init__(self, smec_start_pos, smec_velocity,
          nsamples, scan_duration, prescan_duration)

        # read the velocity jitter spectrum
        jitterfile = os.path.join(os.path.dirname(__file__), 'data',
          jitterspectrum)
        hdulist = pyfits.open(jitterfile)
        hdulistData = hdulist[1].data

        # freq in Hz
        self.freq = hdulistData.frequency
        # psectra in um/s/sqrt(Hz)
        self.specForward = hdulistData.specForward 
        self.specReverse = hdulistData.specReverse 
        self.specMean = hdulistData.specMeanFwRv

        # initialize random number generator
        np.random.seed(seed)

        self.demand_spec = {}

    def _calculate_vel_errors(self):
        # Velocity errors will be derived for scan samples only.
        # Errors for prescan samples will be 0.

        # calculate the frequency axis of the error spectrum
        # that corresponds to the valid scan samples.
        times = self.scan_relative_time[self.prescan_nsamples:]

        demand_df = 1.0 / (times[-1] - times[0])
        if len(times)%2 != 0:
            raise Exception, 'expecting time sequence with even number of points'

        demand_nf = int(len(times) / 2)

        if (demand_df, demand_nf) not in self.demand_spec.keys():
#            print 'creating spectrum for', (demand_df, demand_nf)
            demand_f = np.arange(demand_nf) * demand_df
            demand_spec = np.zeros(np.shape(demand_f))

            for i,f in enumerate(demand_f):
                f_index = (f - self.freq[0]) / (self.freq[1] - self.freq[0])
                ibelow = int(f_index)
                iabove = ibelow + 1
                demand_spec[i] = self.specMean[ibelow] + \
                  (self.specMean[iabove] - self.specMean[ibelow]) * \
                  (f_index - float(ibelow))
    
            self.demand_spec[(demand_df, demand_nf)] = (
              demand_f, demand_spec)
        else:
#            print 'reusing spectrum for', (demand_df, demand_nf)
            demand_f, demand_spec = self.demand_spec[(demand_df, demand_nf)]

        deltafreq = demand_f[1] - demand_f[0]
        max_time = 1.0 / deltafreq
        nfreq = len(demand_f)

        # create random phase using the profile
#        randomphase = np.random.uniform(low=-np.pi, high=np.pi,
#          size=nfreq)    
        randomphase = np.arange(nfreq) * np.pi / float(nfreq)    

        # multiply spectrum by random phase
        mvspec = demand_spec * (np.cos(randomphase) + 
          1j * np.sin(randomphase))

        # convert spectrum to m/s
#        mvspec *= (np.sqrt(deltafreq) * 1.0e-6)

        # butterfly mvspec (real part even, complex part odd)
        btf_mvspec = np.zeros([2*nfreq], np.complex)
        btf_mvspec[:nfreq] = mvspec
        btf_mvspec[nfreq+1:] = np.conj(mvspec)[nfreq-1:0:-1]
        # renormalise as now doubled in length 
        btf_mvspec *= 0.5

        # vel error array covers prescan samples plus scan samples
        vel_error = np.zeros(np.shape(self.scan_relative_time))

        # inverse fft to transform mvspec to time domain. Don't want
        # normalisation so multiply result by 2*nfreq, the length of
        # btf_mvspec. Place result into the valid sample part of the 
        # vel error array
        vel_error[self.prescan_nsamples:] = \
          float(2*nfreq) * np.fft.ifft(btf_mvspec).real

#        plt.figure()
#        plt.plot(errorRA.real)
#        print errorRA.real[:10]
#        plt.plot(errorRA.imag, color='r')
#        plt.savefig('real.png')
#        plt.close()
#        x=1/0

#        return vel_error
        return btf_mvspec.real

    def _calculate_vel_errors_2(self):
        # Velocity errors will be derived for scan samples only.
        # Errors for prescan samples will be 0.

        # This version use corkscrew phase slowly varying over
        # power spectrum. This gives increased errors at scan ends
        # but same result might be obtained by fixing




        # phase across resonances as Gibion recommends.

        deltafreq = self.freq[1] - self.freq[0]
        max_time = 1.0 / deltafreq
        nfreq = len(self.freq)

        # create random phase using the profile
        randomphase = (np.random.random() * np.pi) + (
          np.arange(nfreq) * 3.0 * np.pi / float(nfreq))    

        # multiply spectrum by random phase
        mvspec = self.specMean * (np.cos(randomphase) + 
          1j * np.sin(randomphase))

        # convert spectrum to m/s
        mvspec *= np.sqrt(deltafreq) * 1.0e-6

        # butterfly mvspec (real part even, complex part odd)
        btf_mvspec = np.zeros([2*nfreq], np.complex)
        btf_mvspec[:nfreq] = mvspec
        btf_mvspec[nfreq+1:] = np.conj(mvspec)[nfreq-1:0:-1]
        # renormalise as now doubled in length 
        btf_mvspec *= 0.5

        # inverse fft to transform mvspec to time domain. Don't want
        # normalisation so multiply result by 2*nfreq, the length of
        # btf_mvspec.

        vel_error = float(2*nfreq) * np.fft.ifft(btf_mvspec).real
        # fudge factor to get velocities to roughly match those
        # in Gibion's powerpoint slide
        vel_error *= 7.0

#        vel_error[:5] = 0.0
#        vel_error[-5:] = 0.0

        interp_vel_error = np.zeros([self.nsamples])

        # vel error array covers prescan samples plus scan samples
        scan_vel_error = np.zeros([self.total_nsamples])

        for i in range(self.nsamples):
            i_index = (float(i) / float(self.nsamples)) * float(len(vel_error))
            ibelow = int(i_index)
            iabove = ibelow + 1
            interp_vel_error[i] = vel_error[ibelow] + \
              (vel_error[iabove] - vel_error[ibelow]) * \
              (i_index - float(ibelow))

        # interpolate into valid points of scan array
        scan_vel_error[self.prescan_nsamples:] = interp_vel_error

#        plt.figure()
#        plt.plot(scan_vel_error)
#        plt.savefig('real.png')
#        plt.close()
#        x=1/0

        return scan_vel_error

    def run(self, nscans):
        """Calculate an array with actual positions for the given time
        sequence.
        """
        if nscans%2 != 0:
            raise Exception, 'nscans not multiple of 2'

        times = np.zeros([nscans * self.total_nsamples])
        demand_smec_position = np.zeros([nscans * self.total_nsamples])
        actual_smec_position = np.zeros([nscans * self.total_nsamples])
        smec_integrate = np.zeros([nscans * self.total_nsamples], np.bool)

        # first scan will be forward
        direction = 1.0
        scan_start_time = 0.0
        smec_start_pos = self.smec_start_pos

        for scan in range(nscans):
            # calculate velocity errors for this scan
            vel_error = self._calculate_vel_errors_2()

            self._calculate_smec_position(scan, scan_start_time, smec_start_pos, 
              direction, times, demand_smec_position, actual_smec_position,
              smec_integrate, vel_error)

            # swap direction for next scan
            direction = -direction
            scan_start_time = times[(scan+1)*self.total_nsamples-1]
            smec_start_pos = demand_smec_position[(scan+1)*self.total_nsamples-1]

        return times, demand_smec_position, actual_smec_position, smec_integrate

    def __repr__(self):
        return 'smecposition.HerschelErrors'

