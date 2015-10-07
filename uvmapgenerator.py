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

        if pattern.lower() == 'fixed arc length':
            # the uv pattern is a set of rings, with points distributed
            # along the ring arcs at roughly equal spacing for all rings

            bmax = interferometer['bmax [m]'][row]
            bmin = interferometer['bmin [m]'][row]
            bstep = interferometer['bstep [m]'][row]
            bmin_ang = interferometer['bminAngle[deg]'][row]
            # the time to be spent at each baseline position
            baseline_dwell = interferometer['period [s]'][row]

            # number of baseline radii
            if bmax == bmin:
                n_radii = 1
            else:
                n_radii = int((bmax - bmin) / bstep) + 1

            # range of baseline radii
            b = np.linspace(bmin, bmax, n_radii)

            # length of arc corresponding to arc angle at shortest radius.
            # This will be the arc length used for all rings    
            arcl = (bmin / 2) * np.deg2rad(bmin_ang)

            # get angle to which arclength corresponds for each ring,
            # then number of baselines on ring
            bang_step = np.ones([n_radii])
            n_baselines_rad = np.ones([n_radii], dtype=np.int)

            for i in range(n_radii):
                bang_step[i] = (2 * arcl) / b[i]
                n_baselines_rad[i] = int(((2 * np.pi / bang_step[i]) + 1) / 2)
    
            n_baselines = np.sum(n_baselines_rad)

            # construct the baseline positions
            bxby = collections.OrderedDict()
            t = 0.0

            for i in range(n_radii):
                for j in range(n_baselines_rad[i]):

                    # baseline positions
                    bx = b[i] * np.cos((j-1) * bang_step[i])
                    by = b[i] * np.sin((j-1) * bang_step[i])

                    # timeline for each position
                    # first, set the end of the previous 'inter baseline'
                    # period
                    flag = True
                    bxby[t] = (bx, by, flag)
                    t += 0.001

                    # start and end times at this baseline
                    flag = False
                    bxby[t] = (bx, by, flag)
                    t += baseline_dwell
                    bxby[t] = (bx, by, flag)
                    t += 0.001

                    # start the next 'inter baseline period'
                    flag = True
                    bxby[t] = (bx, by, flag)

                    # 60 seconds between baseline positions
                    t += 60.0

            self.result['pattern'] = pattern
            self.result['bxby'] = bxby
            self.result['bmin'] = bmin
            self.result['bmax'] = bmax

        elif pattern.lower() == 'const l spiral':
            # spiral that conserves angular momentum as the baseline is 
            # increased. If moving from bmin=10m to bmax=100m then 
            # there is a factor of 100 change in angular velocity.

            bmax = interferometer['bmax [m]'][row]
            bmin = interferometer['bmin [m]'][row]
            bstep = interferometer['bstep [m]'][row]
            bmax_period = interferometer['period [s]'][row]

            # baseline increases by bstep for each circuit of spiral
            n_laps = (bmax - bmin) / bstep

            # angular velocity and momentum at bmax 
            omega_start = 2.0 * np.pi / bmax_period
            ang_mom = omega_start * (bmax / 2.0)**2

            bxby = collections.OrderedDict()

            t = 0.0
            b = bmax
            pa = 0.0
            omega = omega_start
            looping = True
            
            while looping:
                bx = b * np.sin(pa)
                by = b * np.cos(pa) 
                flag = False

                bxby[t] = (bx, by, flag)

                # following is probably not the most accurate way to
                # increment the variables, but should be good enough
                t += 1.0
                pa += omega
                b -= (bstep / bmax_period)
                omega = ang_mom / (b / 2.0)**2  

                looping = b > bmin

            self.result['pattern'] = pattern
            self.result['bxby'] = bxby
            self.result['bmin'] = bmin
            self.result['bmax'] = bmax

        elif pattern.lower() == 'spiral':
            bmax = interferometer['bmax [m]'][row]
            bmin = interferometer['bmin [m]'][row]
            n_baselines = int(interferometer['Num Baselines'][row])
            baseline_dwell = interferometer['period [s]'][row]
            bstep = interferometer['bstep [m]'][row]

            # baseline increases by bstep for each circuit of spiral
            n_laps = (bmax - bmin) / bstep

            bxby = collections.OrderedDict()
            t = 0.0

            # n_baselines is total number of points along spiral
            for ib in range(n_baselines):
                # phi is angle around spiral, r its radius
                phi = n_laps * 2 * np.pi * (n_baselines-ib) / n_baselines
                r = bmin + bstep * phi / (2 * np.pi)

                bx = r * np.cos(phi)
                by = r * np.sin(phi)

                # first, set the end of the previous 'inter baseline'
                # period
                flag = True
                bxby[t] = (bx, by, flag)
                t += 0.001

                # start and end times at this baseline
                flag = False
                bxby[t] = (bx, by, flag)
                t += baseline_dwell
                bxby[t] = (bx, by, flag)
                t += 0.001

                # start the next 'inter baseline period'
                flag = True
                bxby[t] = (bx, by, flag)

                # 60 seconds between baseline positions
                t += 60.0

            self.result['pattern'] = pattern
            self.result['bxby'] = bxby
            self.result['bmin'] = bmin
            self.result['bmax'] = bmax

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

            bxby[:, 0] = x1_0 * 2
            bxby[:, 1] = y1_0 * 2

            self.result['pattern'] = pattern
            self.result['bxby'] = bxby
            self.result['bmin'] = R2-R1
            self.result['bmax'] = R1+R2
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


