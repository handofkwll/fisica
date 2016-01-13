"""This module contains classes and methods used to calculate the
observation timeline.
"""

from __future__ import absolute_import

import collections
import numpy as np

import pointing
import smecposition

# This namedtuple will hold the instrument configuration at each
# timestamp in the observation.
Config = collections.namedtuple('Config', [
    'scan_number',
    'time',
    'baseline_x', 'baseline_y', 'baseline_z', 'baseline_flag',
    'smec_position', 'smec_nominal_position',
    'flag', 'smec_vel_error', 
    'pointing1_x', 'pointing1_y',
    'pointing2_x', 'pointing2_y',
    'data', 'pure_data', 'cr_data', 'detector_noise_data'],
    verbose=False)

def find_unflagged_sections(flag):
    """Routine to find ranges [start:end] of sections
    in the given array that are set False. Range means
    that points from flag[start] up to but not including
    flag[end] are False.

    Keyword parameters:
    flag - 1-d boolean array
 
    Returns:
    A list of index ranges [start:end] describing where flag==False
    """

    unflagged_sections = []
    unflagged_section = [-1,-1]
    in_unflagged_section = False
    for i,f in enumerate(flag):
        if not f:
            if in_unflagged_section:
                pass
            else:
                unflagged_section[0] = i
                in_unflagged_section = True
        else:
            if in_unflagged_section:
                unflagged_section[1] = i
                unflagged_sections.append(unflagged_section)
                unflagged_section = [-1,-1]
                in_unflagged_section = False
            else:
                pass

    # tidy up ending
    if in_unflagged_section:
        unflagged_section[1] = i + 1
        unflagged_sections.append(unflagged_section)

    return unflagged_sections


class TimeLineGenerator(object):
    """Class to generate timeline of the simulated observation.

    Contains methods:
    __init__
    run
    __repr__
    """

    def __init__(self, previous_results):
        """Constructor.
        
        Keyword parameters:
            previous_results - dict containing ongoing results structure
                               for the simulation.
        """
        self.previous_results = previous_results
        self.result = collections.OrderedDict()

    def run(self):
        # print 'TimeLineGenerator.run'

        # The FTS mirror positions per scan. This is found from the 
        # FTS results that have already been calculated.
        fts = self.previous_results['fts']
        self.result['smec_start'] = smec_start = fts['smec_start'] / 100.0
        self.result['smec_velocity'] = smec_velocity = fts['smec_velocity'] / 100.0
        self.result['fts_nsample'] = smec_nsample = fts['fts_nsample']
        self.result['smec_scan_duration'] = smec_scan_duration = \
          fts['smec_scan_duration']
        self.result['smec_interscan_duration'] = smec_interscan_duration = \
          fts['smec_interscan_duration']
        self.result['smec_verror_type'] = smec_verror_type = \
          fts['smec_verror_type']
        self.result['nscans'] = nscans = fts['nscans']
        
        # the baselines observed, from 'uvmapgenerator'
        uvmap = self.previous_results['uvmapgenerator']
        bxby = uvmap['bxby']

        # collector parameters
        telescope = self.previous_results['telescope']
        self.result['c1_pointing_error_type'] = \
          c1_pointing_error_type = telescope['c1_pointing_error_type']
        self.result['c2_pointing_error_type'] = \
          c2_pointing_error_type = telescope['c2_pointing_error_type']

        # Construct the obs timeline.
        # The timeline will hold the baseline and FTS mirror position
        # for each time that the detector is read. 

        # The observation consists of a series of FTS scans. The 
        # instrument baseline, pointing etc. can change continuously
        # as the scan is performed.
 
        # Basic timeline is:
        #    Spacecraft baseline sequence gives the length of the 
        #    timeline.
        #
        #    Each baseline contains n_scans FTS forward/back scans,
        #    separated by an 'interscan' period to allow for
        #    acceleration/deceleration etc.

        obs_timeline = {}
        baseline_start_time = 0.0
        inter_baseline_time = 3600.0

        # objects to use for generating pointing errors.
        if 'HERSCHEL' in c1_pointing_error_type.upper():
            pointing1 = pointing.HerschelErrors()
        elif 'ZERO' in c1_pointing_error_type.upper():
            pointing1 = pointing.ZeroErrors()
        else:
            print 'c1 pointing error mode not understood, defaulting to Zero errors'
            pointing1 = pointing.ZeroErrors()

        if 'HERSCHEL' in c2_pointing_error_type.upper():
            pointing2 = pointing.HerschelErrors()
        elif 'ZERO' in c2_pointing_error_type.upper():
            pointing2 = pointing.ZeroErrors()
        else:
            print 'c2 pointing error mode not understood, defaulting to Zero errors'
            pointing2 = pointing.ZeroErrors()

        # object to use for generating SMEC positions
        if 'HERSCHEL' in smec_verror_type.upper():
            smec = smecposition.HerschelErrors(smec_start, smec_velocity,
              smec_nsample, smec_scan_duration, smec_interscan_duration)   
        elif 'ZERO' in smec_verror_type.upper():
            smec = smecposition.ZeroErrors(smec_start, smec_velocity,
              smec_nsample, smec_scan_duration, smec_interscan_duration)   
        else:
            print 'SMEC error mode not understood, defaulting to Zero errors'
            smec = smecposition.ZeroErrors(smec_start, smec_velocity,
              smec_nsample, smec_scan_duration, smec_interscan_duration)   

        bxby_flag = np.array([v[2] for v in bxby.values()])

        # get ranges of times covering unflagged stretches of baseline 
        # position
        unflagged_sections = find_unflagged_sections(bxby_flag)

        # iterate through unflagged stretches

        start_scan = 0

        for section in unflagged_sections:
            section_times = bxby.keys()[section[0]:section[1]]
            start_time = section_times[0]
            end_time = section_times[-1]

            # if there are previous measurements then fill in the times
            # between the last timestamp and the start of this section
            # with null measurements. This gives a roughly continuous
            # timeline that it's easier to add noise and the detector
            # time response to.
            if obs_timeline:
                obs_times = obs_timeline.keys()
                obs_times.sort()
                dt = obs_times[-1] - obs_times[-2]

                t = obs_times[-1] + dt
                while t < start_time:
                    config = Config(-1, t,
                      0.0, 0.0, 0.0, True,
                      obs_timeline[obs_times[-1]].smec_position,
                      obs_timeline[obs_times[-1]].smec_nominal_position,
                      True, 0.0, 
                      0.0, 0.0, 0.0, 0.0, 
                      0.0, 0.0, 0.0, 0.0)

                    obs_timeline[t] = config
                    t += dt

            # calculate number of scans if FTS is run continuously,
            # force to be even
            nscans = (end_time - start_time) / \
              (smec_scan_duration + smec_interscan_duration)
            nscans = int(nscans)
            if nscans%2 == 1:
                nscans -= 1

            # generate actual positions
            smec_time, smec_nominal_position, smec_position,\
              smec_flag, smec_vel_error, scan_number = smec.run(
              nscans, start_scan)
            pointing1_x, pointing1_y = pointing1.run(times=smec_time)
            pointing2_x, pointing2_y = pointing2.run(times=smec_time)

            smec_time += start_time

            itime = 0
            nextrapolated = 0
            for i,t in enumerate(smec_time):

                # interpolate baseline for each smec position
                if t >= section_times[itime+1]:
                    # watch for extrapolation due to rounding errors near
                    # end of time sequence
                    if itime + 2 < len(section_times):
                        itime += 1
                    else:
                        nextrapolated += 1

                lo_time = section_times[itime]
                hi_time = section_times[itime+1]
                bx = bxby[lo_time][0] + \
                  (bxby[hi_time][0] - bxby[lo_time][0]) * \
                  (t - lo_time) / (hi_time - lo_time)       
                by = bxby[lo_time][1] + \
                  (bxby[hi_time][1] - bxby[lo_time][1]) * \
                  (t - lo_time) / (hi_time - lo_time)       
                bflag = bxby[lo_time][2] or bxby[hi_time][2]       

                ib = 0
                config = Config(scan_number[i], t, bx, by, 0.0, False,
                  smec_position[i], smec_nominal_position[i],
                  smec_flag[i], smec_vel_error[i],
                  pointing1_x[i], pointing1_y[i],
                  pointing2_x[i], pointing2_y[i],
                  0.0, 0.0, 0.0, 0.0)

                # sanity check for duplicate times, this should
                # not happen
                if obs_timeline.has_key(t):
                    print 'duplicate time', i, t, ib

                obs_timeline[t] = config

            # increment scan number ready for next group of scans
            start_scan = np.max(scan_number) + 1

            # print warning if any extrapolating was done. This is caused
            # by rounding errors so the number should be 0 or 1, anything
            # larger and you should check that something else is not
            # awry
            if nextrapolated > 0:
                print '..number of baseline points extrapolated %s' % \
                  nextrapolated
        self.result['obs_timeline'] = obs_timeline

        return self.result

    def __repr__(self):
        return '''
TimeLineGenerator:
  C1 pointing errors  : {c1_pointing_errors}
  C2 pointing errors  : {c2_pointing_errors}
  FTS samples/scan    : {fts_nsample}
  timelength length   : {timeline_len}
'''.format(
          c1_pointing_errors=self.result['c1_pointing_error_type'],
          c2_pointing_errors=self.result['c2_pointing_error_type'],
          fts_nsample=self.result['fts_nsample'],
          timeline_len=len(self.result['obs_timeline']))


