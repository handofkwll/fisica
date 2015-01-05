import os.path

import numpy as np
from scipy.io.idl import readsav


class FtsSmec(object):
    """FTS mirror mechanism position simulator.
    """

    def __init__(self, seed=57, jitterdata='SpireVelocityJitter.sav'):
        print 'file', __file__
        self.speed = 498.6982      # cm/s
        self.scale = 1.0

        jitterfile = os.path.join(os.path.dirname(__file__), 'data', jitterdata)
        print 'jitterfile', jitterfile
        data = readsav(jitterfile, verbose=False)
        print 'data', data

        self.freq = data['mydataarray']
        print 'freq', self.freq
        self.freq = data['mydataarray'].freq[0]
        self.spec = data['mydataarray'].spec[0]
        print 'freq', self.freq
        print 'spec', self.spec

        self.deltafreq = self.freq[1] - self.freq[0]

        # initialize random number generator
        np.random.seed(seed)

        self.mpd_start = {1:0.0, -1:4.0}
        self.direction = 1
       

    def actualpos(self):
        """Calculate an array with actual positions for one mirror scan.
        """

        # create random phase using the profile
        randphase = np.random.uniform(low=-np.pi, high=np.pi,
          size=len(self.spec))    

        # multiply spec by random phase
        mvspec = self.spec * (np.cos(randphase) + 1j * np.sin(randphase))

        # normalise
        mvspec *= np.sqrt(self.deltafreq)

        # butterfly mvspec (real part even, complex part odd)
        # ..'centre' of output array is last point in mvspec
        m = len(mvspec)
        btf_mvspec = np.zeros([2*m-1], np.complex)
        btf_mvspec[:m] = mvspec
        btf_mvspec[m:] = np.conj(mvspec)[m-2::-1]

        # inverse fft to transform mvspec to time domain. Don't want
        # normalisation so multiply result by m.
        verror = float(m) * np.fft.ifft(btf_mvspec).real
        print 'verror', verror

        # convert from microns to metres.
        verror *= 1.0e-6

        nv = len(verror)
        velocity = np.ones([nv]) * self.speed * self.direction
        actualvel = velocity + verror

        nyquist = 1.0 / self.deltafreq
        timearray = np.arange(len(verror)) * nyquist / len(verror)
        dt = timearray[1] - timearray[0]

        # compute the mirror path difference
        mpd = np.zeros([nv])
        mpd[0] = self.mpd_start[self.direction]
        for i in range(nv)[1:]:
            mpd[i] = mpd[i-1] + actualvel[i-1] * dt
              
        # reverse direction for next scan
        self.direction *= -1        
