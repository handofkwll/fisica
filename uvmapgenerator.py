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

            bxby = np.zeros([n_baselines,2])     

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

        elif pattern.lower() == 'spiro':

            bmax = interferometer['bmax [m]'][row]
            bmin = interferometer['bmin [m]'][row]
            n_baselines = int(interferometer['Num Baselines'][row])
            bxby = np.zeros([n_baselines, 2])
            R1 = 18
            R2 = 22

            om_tru = 10e-4
            om_tel = om_tru * 10

            nsteps = n_baselines
            dt = 60

            t1_x = R2 * np.cos(np.linspace(0, nsteps-1, nsteps) * dt * om_tru)
            t1_y = R2 * np.sin(np.linspace(0, nsteps-1, nsteps) * dt * om_tru)

            t2_x = -R2 * np.cos(np.linspace(0, nsteps-1, nsteps) * dt * om_tru)
            t2_y = -R2 * np.sin(np.linspace(0, nsteps-1, nsteps) * dt * om_tru)

            x1_0 = t1_x + R1 * np.cos(np.linspace(0, nsteps-1, nsteps) * om_tel * dt)
            y1_0 = t1_y + R1 * np.sin(np.linspace(0, nsteps-1, nsteps) * om_tel * dt)
            x2_0 = t2_x + R1 * np.cos(np.linspace(0, nsteps-1, nsteps) * om_tel * dt)
            y2_0 = t2_y + R1 * np.sin(np.linspace(0, nsteps-1, nsteps) * om_tel * dt)

            #print x1_0
            bxby[:, 0] = x1_0 * 2
            bxby[:, 1] = y1_0 * 2

            self.result['pattern'] = pattern
            self.result['bxby'] = bxby
            self.result['bmin'] = R2-R1
            self.result['bmax'] = R1+R2
            self.result['n_baselines'] = n_baselines
            #print self.result
        else:
            raise Exception, 'unknown baseline pattern: %s' % pattern

        return self.result

    def __repr__(self):
        return '''
UVMapGenerator:
  uv pattern      : {pattern}
'''.format(
          pattern=self.result['pattern'])


