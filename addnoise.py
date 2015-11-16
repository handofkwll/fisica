from __future__ import absolute_import

import collections
import numpy as np

def addDetectorTimeConstant(obs_timeline, time_constant):
    """This function convolves the data stream with exp(-t/tau) (t>0)
    where tau is the detector time constant.
    """
    if time_constant < 0.0:
        raise Exception, 'invalid detector time constant: %s' % time_constant
    elif time_constant == 0.0:
        pass
    else:
        # get timestamps and data from observation timeline
        obs_times = obs_timeline.keys()
        obs_times.sort()
        obs_times = np.array(obs_times)

        # construct data timeline and detector response
        data = []
        for t in obs_times:
            data.append(obs_timeline[t].data)
        data = np.array(data)
        detector = np.exp(-(obs_times - obs_times[0]) / time_constant)

        # uncomment the following to put spike into time series and test detector
        # constant
        #data[:] = 0.0
        #data[10] = 1.0

        # apply the detector time constant through Fourier convolution
        ftdata = np.fft.fft(data)
        ftdetector = np.fft.fft(detector)
        ftdata *= ftdetector
        data = np.fft.ifft(ftdata)
        data = data.real

        # set the modified values
        for it,t in enumerate(obs_times):
            config = obs_timeline[t]
            config = config._replace(data = data[it])
            obs_timeline[t] = config

def addGlitches(obs_timeline, timesOfGlitches, widthsOfGlitches,
  peaksOfGlitches):
    """This function generates the profiles for the CR glitches and puts them 
    at predetermined positions on the detector timeline.
    """
    obs_times = obs_timeline.keys()
    obs_times.sort()
    obs_times = np.array(obs_times)

    for nglitch,glitch_time in enumerate(timesOfGlitches):
        glitchSpread = obs_times[(obs_times > (glitch_time - 0.02)) & 
          (obs_times < (glitch_time + 0.06))]

        for t in glitchSpread:
            peak = peaksOfGlitches[nglitch]
            location = glitch_time
            width = widthsOfGlitches[nglitch]

            # Moyal function
            glitchdata = peak * np.exp(-0.5 * (
              ((t - location) / width) + np.exp(-(t - location) / width)
              )) / np.sqrt(2.0 * np.pi)

            data = obs_timeline[t].data
            data += glitchdata

            # set the measured value
            config = obs_timeline[t]
            config = config._replace(data = data, cr_data=glitchdata)
            obs_timeline[t] = config

def addNoise(obs_timeline, NEPtot, NEPsky, NEPdet, efficiency,
  f_acq, knee_freq):
    """This function generates random noise corresponding to the 
    total NEP, folds in the effect of 1/f noise and adds the
    result to the detector timeline.
    """
    obs_times = obs_timeline.keys()
    obs_times.sort()
    obs_times = np.array(obs_times)

    # construct data timeline
    data = []
    for t in obs_times:
        data.append(obs_timeline[t].data)
    data = np.array(data)

    # get overall NEP
    NEPall = np.sqrt(NEPtot**2 + NEPsky**2 + NEPdet**2) / efficiency
    sigma = NEPall * np.sqrt(f_acq)
    noise = sigma * np.random.standard_normal(data.shape)

    # 1/f noise
    # ..next line works if obs_times is a regular arithemtic progression
    #   (which it should be, can uncomment next 2 lines to check)
    #print np.min(obs_times[1:] - obs_times[:-1])
    #print np.max(obs_times[1:] - obs_times[:-1])
    freq_vec = np.fft.fftfreq(len(obs_times), obs_times[1] - obs_times[0])
    # ..work around fact that freq_vec[0] = 0 
    freq_vec[0] = freq_vec[1]
    f_noise = (1 + (knee_freq / np.abs(freq_vec)))**0.5
     
    noise = np.fft.ifft(np.fft.fft(noise) * f_noise)
    noise = noise.real

    # set the modified values
    for it,t in enumerate(obs_times):
        config = obs_timeline[t]
        config = config._replace(data = data[it] + noise[it])
        obs_timeline[t] = config

def generatePeaks():
    """This function generates random numbers based on the Moyal distribution. 
    The 'Rejection Method' from http://nbviewer.ipython.org/github/unpingco/
    Python-for-Signal-Processing/blob/master/Sampling_Monte_Carlo.ipynb
    is used.
    The code is adapted from an original version written by G.Makiwa, 
    U.Lethbridge, Canada.
    """
    # p contains the parameters describing the Moyal distribution
    # of the cosmic ray amplitudes from Spire FTS data [peak, offset, width]
    p =   [1.69053134e+02, 4.78682617e-05, 3.65722650e-05]

    # f is the Moyal distribution (no normalised)
    f = lambda x: np.exp(-0.5 * (((x - p[1]) / p[2]) + np.exp(-(x - p[1]) / p[2]))) /\
      (np.sqrt(2.0 * np.pi))
    
    # u1 is a list of 10000 random numbers scaled between 0 and the max x
    # of the sample range
    u1 = np.random.rand(10000) * 0.002
    # M normalises the Moyal distribution to run between 0 and 1,
    # u2 is 10000 random samples between 0 and 1,
    # idx are the indices of points in u2 where those > corresponding
    # f(u1)/M are rejected,
    # myRandomPeaks are the sample points indexed
    M = 0.242
    u2 = np.random.rand(10000)
    idx=np.where(u2<=f(u1)/M)[0] # rejection criterion
    myRandomPeaks = u1[idx]

    # the Moyal distribution itself
    x = np.linspace(0.000, 0.002, 10000)
    moyal = f(x)
    
    return myRandomPeaks, moyal, x

def generateRandomGlitches(obs_timeline, cr_rate):
    """Given the length of the observed timeline, this function 
    returns the following information:
    - numberOfGlitches: number of cosmic ray (CR) hits in the detector
                        timeline. On average there are 2.5 CR hits 
                        in 66.6s.
    - timesOfGlitches : sample times when CR hits occur 
    - widthsOfGlitches: widths of CR glitches 
    - peaksOfGlitches: peaks of CR glitches                 
    """
    obs_times = obs_timeline.keys()
    obs_times.sort()
    timespan = obs_times[-1] - obs_times[0]

    numberOfGlitches = np.ceil(timespan * cr_rate)
    timesOfGlitches = np.random.uniform(0, timespan, numberOfGlitches)
    myRandomPeaks, moyal, x = generatePeaks()
    peaksOfGlitches = myRandomPeaks[:numberOfGlitches]

    # All Spire FTS glitches had similar widths: 0.00097
    widthsOfGlitches = np.ones(numberOfGlitches) * 0.00097
    
    #print 'timespan =', timespan
    #print 'cr rate =', cr_rate
    #print 'numberOfGlitches =', numberOfGlitches 
    #print 'timesOfGlitches =', timesOfGlitches
    #print 'widthsOfGlitches =', widthsOfGlitches
    #print 'peaksOfGlitches  =', peaksOfGlitches

    return numberOfGlitches, timesOfGlitches, widthsOfGlitches, peaksOfGlitches,\
      moyal, x

def nepSky(obs_timeline, cr_rate):
    """Given the length of the observed timeline, this function 
    returns the following information:
    - numberOfGlitches: number of cosmic ray (CR) hits in the detector
                        timeline. On average there are 2.5 CR hits 
                        in 66.6s.
    - timesOfGlitches : sample times when CR hits occur 
    - widthsOfGlitches: widths of CR glitches 
    - peaksOfGlitches: peaks of CR glitches                 
    """
    obs_times = obs_timeline.keys()
    obs_times.sort()
    timespan = obs_times[-1] - obs_times[0]

    numberOfGlitches = np.ceil(timespan * cr_rate)
    timesOfGlitches = np.random.uniform(0, timespan, numberOfGlitches)
    myRandomPeaks, moyal, x = generatePeaks()
    peaksOfGlitches = myRandomPeaks[:numberOfGlitches]

    # All Spire FTS glitches had similar widths: 0.00097
    widthsOfGlitches = np.ones(numberOfGlitches) * 0.00097
    
    #print 'timespan =', timespan
    #print 'cr rate =', cr_rate
    #print 'numberOfGlitches =', numberOfGlitches 
    #print 'timesOfGlitches =', timesOfGlitches
    #print 'widthsOfGlitches =', widthsOfGlitches
    #print 'peaksOfGlitches  =', peaksOfGlitches

    return numberOfGlitches, timesOfGlitches, widthsOfGlitches, peaksOfGlitches,\
      moyal, x


class AddNoise(object):
    """Class to add simulated noise, cosmic rays and detector
       time response to the simulated observation.
    """

    def __init__(self, parameters, previous_results):
        """Constructor.
        
        Parameters:
            previous_results - dict containing ongoing results structure
                               for the simulation.
        """
        self.parameters = parameters
        self.previous_results = previous_results
        self.result = collections.OrderedDict()

    def run(self):
        #print 'AddNoise.run'

        # get relevant parameters read from Excel file
        detectors = self.previous_results['loadparameters']['substages']\
          ['Detectors']
        row = detectors['Detector time constant'].keys()[0]
        detector_time_constant = detectors['Detector time constant'][row]
        # convert to seconds
        detector_time_constant *= 1.0e-6

        row = detectors['knee freq'].keys()[0]
        knee_freq = detectors['knee freq'][row]

        row = detectors['Acquisiton frequency (n tau)'].keys()[0]
        f_acq = detectors['Acquisiton frequency (n tau)'][row]

        row = detectors['CosmicRays'].keys()[0]
        cosmic_ray_rate = detectors['CosmicRays'][row]

        backgroundnoise = self.previous_results['backgroundnoise']
        NEPtot = backgroundnoise['NEPtot']
        NEPtot = backgroundnoise['NEPtot']

        detectornoise = self.previous_results['detectornoise']
        NEPdet = detectornoise['NEP']
        eta = detectornoise['eta']

        # get the observation timeline constructed earlier
        timeline = self.previous_results['timeline']
        obs_timeline = timeline['obs_timeline']

        # construct object to use for generating cosmic ray strikes
        nglitches, glitchtimes, glitchwidths, glitchpeaks, moyal,\
          moyalx = generateRandomGlitches(obs_timeline, cosmic_ray_rate)

        # add background and detector noise to data
        print 'dummy NEPsky'
        NEPsky = 0.0
        print NEPtot, NEPsky, NEPdet, eta, f_acq, knee_freq
        addNoise(obs_timeline, NEPtot, NEPsky, NEPdet, eta,
          f_acq, knee_freq)

        # add the cosmic ray strikes
        addGlitches(obs_timeline, glitchtimes, glitchwidths,
          glitchpeaks)

        # apply the detector time constant
        addDetectorTimeConstant(obs_timeline, detector_time_constant)   

        self.result['cosmic ray rate'] = cosmic_ray_rate
        self.result['detector time constant'] = detector_time_constant
        self.result['glitch peaks'] = glitchpeaks
        self.result['glitch peaks moyal'] = (moyal, moyalx)
        self.result['glitch times'] = glitchtimes
        self.result['timeline'] = obs_timeline

        return self.result

    def __repr__(self):
        return '''
AddNoise:
  Cosmic Ray Rate        : {cr_rate}
  Detector time constant : {detector_time_constant}
'''.format(
          cr_rate=self.result['cosmic ray rate'],
          detector_time_constant=self.result['detector time constant'])


