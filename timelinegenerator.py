from __future__ import absolute_import

import collections
import numpy as np

import pointing
import smecposition


class TimeLineGenerator(object):
    """Class to generate timeline framework of the simulated observation.
    """

    def __init__(self, previous_results):
        self.previous_results = previous_results
        self.result = collections.OrderedDict()

    def run(self):
        print 'TimeLineGenerator.run'

        # The FTS mirror positions per scan. This is found from the 
        # FTS results that should have already been calculated.
        fts = self.previous_results['fts']
        self.result['smec_start'] = smec_start = fts['smec_start']
        self.result['smec_velocity'] = smec_velocity = fts['smec_velocity']
        self.result['fts_nsample'] = smec_nsample = fts['fts_nsample']
        self.result['smec_scan_duration'] = smec_scan_duration = \
          fts['smec_scan_duration']
        self.result['smec_interscan_duration'] = smec_interscan_duration = \
          fts['smec_interscan_duration']
        
        # the number of baselines observed, from 'uvmapgenerator'
        uvmap = self.previous_results['uvmapgenerator']
        n_baselines = uvmap['n_baselines']
        bxby = uvmap['bxby']

        # Construct the obs framework that will give one FTS
        # scan (forward/back) at each baseline.
     
        # The framework will hold the baseline and FTS mirror position
        # for each time that the detector is read. 

        # The observation consists of a series of FTS scans. The 
        # instrument baseline, pointing etc. can change continuously
        # as the scan is performed.
 
        # Basic framework is:
        #    n_baselines 'baselines',
        #    separated by a time for spacecraft configuration.
        #
        #    Each baseline contains n_scans FTS forward/back scans,
        #    separated by an 'interscan' period to allow for
        #    accleration/deceleration etc.
        #
        #    Each FTS forward or back scan takes fts_nsamples 
        #    lasting (opd_end - opd_start) / (fts_nsamples * fts_speed)  

        obs_framework = {}
        baseline_start_time = 0.0
        nscans = 4
        inter_baseline_time = 3600.0

        Config = collections.namedtuple('Config', [
          'scan_number',
          'time',
          'baseline_x', 'baseline_y', 'baseline_z', 'baseline_number',
          'smec_position', 'smec_nominal_position',
          'integrate', 
          'pointing1_x', 'pointing1_y',
          'pointing2_x', 'pointing2_y',
          'data'],
          verbose=False)

        # objects to use for generating pointing errors
        pointing1 = pointing.HerschelErrors()
        pointing2 = pointing.ZeroErrors()

        # object to use for generating SMec positions
        smec = smecposition.HerschelErrors(smec_start, smec_velocity,
          smec_nsample, smec_scan_duration, smec_interscan_duration)   
#        smec = smecposition.ZeroErrors(smec_start, smec_velocity,
#          smec_nsample, smec_scan_duration, smec_interscan_duration)   

        for ib in range(n_baselines):
            baseline_start_time += inter_baseline_time

            # generate actual positions
            smec_time, smec_nominal_position, smec_position,\
              smec_integrate = smec.run(nscans)
            pointing1_x, pointing1_y = pointing1.run(times=smec_time)
            pointing2_x, pointing2_y = pointing2.run(times=smec_time)

            smec_time += baseline_start_time            
            for i,t in enumerate(smec_time):
                config = Config(ib, t, bxby[ib][0], bxby[ib][1], 0.0, ib,
                  smec_position[i], smec_nominal_position[i],
                  smec_integrate[i], pointing1_x[i], pointing1_y[i],
                  pointing2_x[i], pointing2_y[i], None)
                obs_framework[t] = config 

        self.result['obs_framework'] = obs_framework

        return self.result

    def __repr__(self):
        return '''
TimeLineGenerator:
  FTS samples/scan: {fts_nsample}
  obs.framework length: {len_framework}
'''.format(
          fts_nsample=self.result['fts_nsample'],
          len_framework=len(self.result['obs_framework']))


