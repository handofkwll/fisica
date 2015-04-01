from __future__ import absolute_import

import numpy as np
import os.path
import astropy.io.fits as pyfits
import scipy.interpolate

import matplotlib.pyplot as plt

class GaussianErrors(object):
    """Class to fill timeline with Gaussian pointing errors.
    """

    def __init__(self, mu0, sigma0, mu1, sigma1):
        self.mu0 = mu0
        self.sigma0 = sigma0
        self.mu1 = mu1
        self.sigma1 = sigma1

    def fill(self, config_list):
        print 'GaussianPointingErrors.fill'
        x0 = self.mu0 + self.sigma0 * np.random.randn(len(config_list))
        y0 = self.mu0 + self.sigma0 * np.random.randn(len(config_list))
        x1 = self.mu1 + self.sigma1 * np.random.randn(len(config_list))
        y1 = self.mu1 + self.sigma1 * np.random.randn(len(config_list))

        for i,config in enumerate(config_list):
            newconfig = config._replace(point1_x_error=x0[i],
              point1_y_error=y0[i],
              point2_x_error=x1[i],
              point2_y_error=y1[i])
            config_list[i] = newconfig

    def __repr__(self):
        return 'pointing.GaussianErrors'


class HerschelErrors(object):
    """Class to create pointing error sequence derived from
    Herschel pointing error spectrum.
    """

    def __init__(self, seed=57, jitterspectrum='PointingJitter1342270191.fits'):

        # read the pointing jitter power spectrum
        jitterfile = os.path.join(os.path.dirname(__file__), 'data',
          jitterspectrum)
        hdulist = pyfits.open(jitterfile)
        hdulistData = hdulist[1].data

        # freq in Hz
        self.freq = hdulistData.frequency
        # power in arcsec/sec/rt(Hz)
        self.jitterRA = hdulistData.specRA
#        # sanity check by setting spectrum to 1Hz sine wave        
#        self.jitterRA = np.zeros(np.shape(self.jitterRA))
#        one_hertz = int(1.0 / (self.freq[1] - self.freq[0])) 
#        self.jitterRA[one_hertz] = 1.0
        self.jitterDec = hdulistData.specDec

        # initialize random number generator
        np.random.seed(seed)

        self.resampledRA = {}
        self.resampledDec = {}

    def run(self, times):
        """Calculate an array with actual positions for the given time
        sequence. Don't resample power spectrum to match desired
        times as this takes too long - instead just calculate a new
        sequence with random phase and interpolate from the 
        calculated times to the required times. 
        """
        deltafreq = self.freq[1] - self.freq[0]
        max_time = 1.0 / deltafreq
        nfreq = len(self.freq)

        # create random phase using the profile
        randomphaseRA = np.random.uniform(low=-np.pi, high=np.pi,
          size=nfreq)    
        randomphaseDec = np.random.uniform(low=-np.pi, high=np.pi,
          size=nfreq)    

        # multiply spectrum by random phase
        mvspecRA = self.jitterRA * (np.cos(randomphaseRA) + 
          1j * np.sin(randomphaseRA))
        mvspecDec = self.jitterDec * (np.cos(randomphaseDec) + 
          1j * np.sin(randomphaseDec))

        # convert spectrum to arcsec/sec
        mvspecRA *= np.sqrt(deltafreq)
        mvspecDec *= np.sqrt(deltafreq)

        # butterfly mvspecs (real part even, complex part odd)
        btf_mvspecRA = np.zeros([2*nfreq-1], np.complex)
        btf_mvspecRA[:nfreq] = mvspecRA
        btf_mvspecRA[nfreq:] = np.conj(mvspecRA)[nfreq-1:0:-1]
        # renormalise as now doubled in length 
        btf_mvspecRA *= 0.5
        btf_mvspecDec = np.zeros([2*nfreq-1], np.complex)
        btf_mvspecDec[:nfreq] = mvspecDec
        btf_mvspecDec[nfreq:] = np.conj(mvspecDec)[nfreq-1:0:-1]
        btf_mvspecDec *= 0.5

        # inverse fft to transform mvspec to time domain. Don't want
        # normalisation so multiply result by 2*nfreq-1, the length of
        # btf_mvspec. Gibion uses [2*nfreq-1] + 1?
        errorRA = float(2*nfreq-1) * np.fft.ifft(btf_mvspecRA).real
        errorDec = float(2*nfreq-1) * np.fft.ifft(btf_mvspecDec).real
#        plt.figure()
#        plt.plot(errorRA.real)
#        print errorRA.real[:10]
#        plt.plot(errorRA.imag, color='r')
#        plt.savefig('real.png')
#        plt.close()
#        x=1/0

        # construct times to match. There are 2*nfreq-1 points running
        # from time=0 to time=max_time 
        errorTimes = np.arange(2*nfreq-1) * max_time / float(2*nfreq-1)

        # now use linear interpolation to go from errorTimes to
        # the requested times 
        deltaErrorTime = errorTimes[1] - errorTimes[0]
        eRA = np.zeros(np.shape(times))
        eDec = np.zeros(np.shape(times))
        for i,t in enumerate(times):
            error_index = (t - errorTimes[0]) / deltaErrorTime
            ibelow = int(error_index)
            iabove = ibelow + 1
            eRA[i] = errorRA[ibelow] + (errorRA[iabove] - errorRA[ibelow]) * \
              (error_index - float(ibelow))
            eDec[i] = errorDec[ibelow] + (errorDec[iabove] - errorDec[ibelow]) * \
              (error_index - float(ibelow))

        return eRA, eDec

    def __repr__(self):
        return 'pointing.HerschelErrors'


class ZeroErrors(object):
    """Class to fill timeline with 0 pointing errors.
    """

    def run(self, times):
#        print 'pointing.ZeroErrors.run'
        
        pointing_x = np.zeros(np.shape(times))
        pointing_y = np.zeros(np.shape(times))

        return pointing_x, pointing_y

    def __repr__(self):
        return 'pointing.ZeroErrors'

