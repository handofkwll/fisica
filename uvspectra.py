from __future__ import absolute_import

import collections
import math
import numpy as np
import pp

import common.commonobjects as co


class UVspectra(object):
    """Class to compute spectra from each interferogram.
    """

    def __init__(self, previous_results, job_server):
        self.previous_results = previous_results
        self.job_server = job_server

        self.uvspectra = []
        self.result = collections.OrderedDict()      

    def run(self):
        print 'Uvspectra.run'

        # get observation list
        observe = self.previous_results['observe']
        obs_framework = observe['observed_framework']

        # and calculate the result for each configuration
        UVspectrum = collections.namedtuple('uvspectrum', [
          'scan_number',
          'baseline_x', 'baseline_y', 'baseline_z', 'baseline_number',
          'interferogram', 'mirror_position', 'spectrum', 'wavenumber'],
          verbose=False)

        fts_position = []
        fts_data = []

        for config in obs_framework:
            if config.fts_start:
                if fts_data:
                    # reduce previous scan - assume symmetric -tau to +tau, so first
                    # shift 0 freq from centre of scan
                    temp = np.fft.fftshift(fts_data)
                    # spectrum is complex
                    spectrum = np.fft.fft(temp)
                    wavenumber = np.fft.fftfreq(n=spectrum.size,
                      d=2.0 * abs(fts_position[1] - fts_position[0])) / 100.0
                    
                    # save it
                    uvspectrum = UVspectrum(config.scan_number,
                      baseline[0], 
                      baseline[1],
                      0.0,
                      config.baseline_number,
                      np.array(fts_data),
                      np.array(fts_position),
                      spectrum,
                      wavenumber)

                    self.uvspectra.append(uvspectrum)

                # start new scan
                fts_position = [config.fts_nominal_position]
                fts_data = [config.data]

            else:

                # add points to current scan
                fts_position.append(config.fts_nominal_position)
                fts_data.append(config.data)
            baseline = (config.baseline_x, config.baseline_y)

        # deal with last scan
        if fts_data:
            temp = np.fft.fftshift(fts_data)
            spectrum = np.fft.fft(temp)
            wavenumber = np.fft.fftfreq(n=spectrum.size,
              d=2.0 * abs(fts_position[1] - fts_position[0])) / 100.0
                    
            uvspectrum = UVspectrum(config.scan_number,
              baseline[0], 
              baseline[1],
              0.0,
              config.baseline_number,
              np.array(fts_data),
              np.array(fts_position),
              spectrum,
              wavenumber)

            self.uvspectra.append(uvspectrum)

        self.result['uvspectra'] = self.uvspectra
        return self.result
                
    def __repr__(self):
        return '''
UVspectra:
  #uvspectra : {num_uvspectra}
'''.format(
          num_uvspectra=len(self.uvspectra))

