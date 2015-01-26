from __future__ import absolute_import

import numpy as np
import os.path
import pyfits

import matplotlib.pyplot as plt


class ZeroErrors(object):
    """Class to generate SMEC positions with no velocity errors.
    """

    def __init__(self, smec_start_pos, smec_velocity, nsamples,
      scan_duration, prescan_duration):
        """Constructor.

        Parameters:
        smec_start_pos   - start position for this scan [m]
        smec_velocity    - nominal smec velocity [m/s] 
        nsamples         - number of samples in scan
        scan_duration    - duration of scan [s]
        prescan_duration - duration of 'prescan' period, where this is 
                           length of time before the scan starts. This to
                           'catch' transients triggered before scan start
                           that will affect the scan e.g. cosmic ray
                           strikes [s].
        """
        self.smec_start_pos = smec_start_pos
        self.smec_velocity = smec_velocity
        self.nsamples = nsamples
        self.scan_duration = scan_duration
        self.prescan_duration = prescan_duration

        # construct sampling times
        sample_time  = scan_duration / float(nsamples)
        total_time = prescan_duration + scan_duration
        self.total_nsamples = int(np.ceil(total_time / sample_time))
        self.prescan_nsamples = self.total_nsamples - nsamples
        # +1 is there to ensure that end/start of consecutive scans do not
        # overlap in time
        self.scan_relative_time = np.arange(self.total_nsamples, dtype=np.float) + 1.0
        self.scan_relative_time *= sample_time

        # smec demand velocity, 0 for prescan samples
        self.smec_vel = np.ones([self.total_nsamples]) * self.smec_velocity
        self.smec_vel[:self.prescan_nsamples] = 0.0

    def _calculate_smec_position(self, scan, scan_start_time, smec_start_pos, 
      direction, relative_time, vel_error, demand_smec_position,
      actual_smec_position, smec_flag):
        """Method to calculate smec position from scan info.
        """

        # times 'internal' to this scan
        relative_time[scan*self.total_nsamples : (scan+1)*self.total_nsamples] = \
          scan_start_time + self.scan_relative_time

        # prescan position for now just smec_start_pos
        offset = scan * self.total_nsamples
        demand_smec_position[offset] = smec_start_pos
        actual_smec_position[offset] = smec_start_pos

        # iterate through the times integrating the velocity
        for i in range(1, len(self.scan_relative_time)):
            demand_smec_position[offset+i] = demand_smec_position[offset+i-1] + \
              (self.scan_relative_time[i] - self.scan_relative_time[i-1]) * \
              self.smec_vel[i-1] * direction

            actual_smec_position[offset+i] = actual_smec_position[offset+i-1] + \
              (self.scan_relative_time[i] - self.scan_relative_time[i-1]) * \
              (self.smec_vel[i-1] + vel_error[i]) * direction

        smec_flag[offset+self.prescan_nsamples:offset+self.total_nsamples] = \
          False

    def run(self, nscans):
        """Method to calculate smec data for 'nscans' scans.
        """
#        print 'ZeroSmecErrors.run'
        if nscans%2 != 0:
            raise Exception, 'nscans not multiple of 2'

        times = np.zeros([nscans * self.total_nsamples])
        demand_smec_position = np.zeros([nscans * self.total_nsamples])
        actual_smec_position = np.zeros([nscans * self.total_nsamples])
        smec_flag = np.ones([nscans * self.total_nsamples], np.bool)
        # smec vel errors are 0.
        smec_vel_error = np.zeros([nscans * self.total_nsamples])
        scan_vel_error = np.zeros([self.total_nsamples])

        # first scan will be forward
        direction = 1.0
        scan_start_time = 0.0
        smec_start_pos = self.smec_start_pos

        # iterate through scans, calculating each in turn
        for scan in range(nscans):
            self._calculate_smec_position(scan, scan_start_time, smec_start_pos, 
              direction, times, scan_vel_error, demand_smec_position,
              actual_smec_position, smec_flag)

            # swap direction for next scan
            direction = -direction
            scan_start_time = times[(scan+1)*self.total_nsamples-1]
            smec_start_pos = demand_smec_position[(scan+1)*self.total_nsamples-1]

        return times, demand_smec_position, actual_smec_position, smec_flag,\
          smec_vel_error

    def __repr__(self):
        return 'smecposition.ZeroErrors'


class HerschelErrors(ZeroErrors):
    """Class to generate SMEC positions with velocity errors
    derived from Herschel data.
    """

    def __init__(self, smec_start_pos, smec_velocity, nsamples,
      scan_duration, prescan_duration, seed=57,
      jitterspectrum='VelocityJitter1342270191.fits'):
        """Constructor.
        Parameters:
        smec_start_pos   - start position for this scan [m]
        smec_velocity    - nominal smec velocity [m/s] 
        nsamples         - number of samples in scan
        scan_duration    - duration of scan [s]
        prescan_duration - duration of 'prescan' period, where this is 
                           length of time before the scan starts. This to
                           'catch' transients triggered before scan start
                           that will affect the scan e.g. cosmic ray
                           strikes [s].
        seed             - seed for random number generator.
        jitterspectrum   - File containing velocity jitter spectrum
                           from Herschel.
        """

        # Instantiate the base class
        ZeroErrors.__init__(self, smec_start_pos, smec_velocity,
          nsamples, scan_duration, prescan_duration)

        # read the velocity jitter spectrum
        jitterfile = os.path.join(os.path.dirname(__file__), 'data',
          jitterspectrum)
        hdulist = pyfits.open(jitterfile)
        hdulistData = hdulist[1].data

        # freq in Hz
        self.freq = hdulistData.frequency
        # spectra in um/s/sqrt(Hz)
        self.specForward = hdulistData.specForward 
        self.specReverse = hdulistData.specReverse 
        self.specMean = hdulistData.specMeanFwRv

        # initialize random number generator
        np.random.seed(seed)

        self.demand_spec = {}

    def _calculate_vel_errors(self):
        # NOT USED AT PRESENT.
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
        mvspec *= (np.sqrt(deltafreq) * 1.0e-6)

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
        # Calculate velocity errors for a scan.
 
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

        # linear interpolate into valid points of scan array
        interp_vel_error = np.zeros([self.nsamples])
        for i in range(self.nsamples):
            i_index = (float(i) / float(self.nsamples)) * float(len(vel_error))
            ibelow = int(i_index)
            iabove = ibelow + 1
            interp_vel_error[i] = vel_error[ibelow] + \
              (vel_error[iabove] - vel_error[ibelow]) * \
              (i_index - float(ibelow))

        # vel error array covers prescan samples plus scan samples
        scan_vel_error = np.zeros([self.total_nsamples])
        scan_vel_error[self.prescan_nsamples:] = interp_vel_error

        return scan_vel_error

    def run(self, nscans):
        """Method to calculate smec data for 'nscans' scans.
        """
        if nscans%2 != 0:
            raise Exception, 'nscans not multiple of 2'

        times = np.zeros([nscans * self.total_nsamples])
        demand_smec_position = np.zeros([nscans * self.total_nsamples])
        actual_smec_position = np.zeros([nscans * self.total_nsamples])
        smec_flag = np.ones([nscans * self.total_nsamples], np.bool)
        smec_vel_error = np.zeros([nscans * self.total_nsamples])

        # first scan will be forward
        direction = 1.0
        scan_start_time = 0.0
        smec_start_pos = self.smec_start_pos

        for scan in range(nscans):
            # calculate velocity errors for this scan
            scan_vel_error = self._calculate_vel_errors_2()
            smec_vel_error[scan*self.total_nsamples:(scan+1)*self.total_nsamples] = \
              scan_vel_error

            # and now the mirror position
            self._calculate_smec_position(scan, scan_start_time, smec_start_pos, 
              direction, times, scan_vel_error, demand_smec_position,
              actual_smec_position, smec_flag)

            # swap direction for next scan
            direction = -direction
            scan_start_time = times[(scan+1)*self.total_nsamples-1]
            smec_start_pos = demand_smec_position[(scan+1)*self.total_nsamples-1]

        return times, demand_smec_position, actual_smec_position, smec_flag,\
          smec_vel_error

    def __repr__(self):
        return 'smecposition.HerschelErrors'

