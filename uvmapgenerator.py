from __future__ import absolute_import

import collections
import numpy as np


class UVMapGenerator(object):
    """Class to generate the UV map of the simulated observation.
    """

    def __init__(self, parameters, previous_results):
        self.parameters = parameters
        self.previous_results = previous_results
        self.result = collections.OrderedDict()

    def run(self):
        #print 'UVMapGenerator.run'

        interferometer = self.parameters['substages']['Interferometer']

        # slightly awkward getting pattern value as it's keyed by its
        # row number in the spreadsheet
        row = interferometer['Pattern'].keys()[0]
        pattern = interferometer['Pattern'][row]

        if pattern.lower() == 'spiral':
            bmax = interferometer['bmax [m]'][row]
            bmin = interferometer['bmin [m]'][row]
            n_baselines = int(interferometer['Num Baselines'][row])
            bstep = interferometer['bstep [m]'][row]

            # baseline increases by bstep for each circuit of spiral
            n_laps = (bmax - bmin) / bstep

#            n_baselines = 50
            bxby = np.zeros([n_baselines,2])     
#            bxby[0,0] = 10.0
#            bxby[0,1] = 0.0
#            bxby[1,0] = 20.0
#            bxby[1,1] = 0.0
#            bxby[2,0] = 30.0
#            bxby[2,1] = 0.0
#            bxby[3,0] = 40.0
#            bxby[3,1] = 0.0
#            bxby[4,0] = 0.0
#            bxby[4,1] = 10.0
#            bxby[5,0] = 0.0
#            bxby[5,1] = 20.0
#            bxby[6,0] = 0.0
#            bxby[6,1] = 30.0
#            bxby[7,0] = 0.0
#            bxby[7,1] = 40.0
#            bxby[8,0] = 40.0
#            bxby[8,1] = 40.0
#            bxby[9,0] = 20.0
#            bxby[9,1] = 20.0
#            bxby[10,0] = 30.0
#            bxby[10,1] = 30.0

            # n_baselines is total number of points along spiral
            for ib in range(n_baselines):
                # phi is angle around spiral, r its radius
                phi = n_laps * 2 * np.pi * (n_baselines-ib) / n_baselines
                r = bmin + bstep * phi / (2 * np.pi)

                bxby[ib,:] = [r * np.cos(phi), r * np.sin(phi)]

            self.result['pattern'] = pattern
            self.result['bxby'] = bxby
            self.result['bmin'] = bmin
            self.result['bmax'] = bmax
            self.result['n_baselines'] = n_baselines

        else:
            raise Exception, 'unknown baseline pattern: %s' % pattern

        return self.result

    def __repr__(self):
        return '''
UVMapGenerator:
  uv pattern      : {pattern}
'''.format(
          pattern=self.result['pattern'])


