from __future__ import absolute_import

import collections
import numpy as np

import common.commonobjects as co

def find_flagged_sections(flag):
    """Routine to find ranges [start:end] of sections
    in the given array that are set true.

    Parameters:
    flag - 1-d boolean array
 
    Returns:
    A list of index ranges [start:end] describing where flag==True
    """

    flagsections = []
    flagsection = [-1,-1]
    inflagsection = False
    for i,f in enumerate(flag):
        if f:
            if inflagsection:
                pass
            else:
                flagsection[0] = i
                inflagsection = True
        else:
            if inflagsection:
                flagsection[1] = i-1
                flagsections.append(flagsection)
                flagsection = [-1,-1]
                inflagsection = False
            else:
                pass

    # tidy up ending
    if inflagsection:
        flagsection[1] = i
        flagsections.append(flagsection)

    return flagsections


class ReduceInterferogram(object):
    """Class to reduce the interferogram scans taken at each 
    baseline position and convert them to a spectrum through
    the uv cube at that point.
    """

    def __init__(self, previous_results, job_server):
        """Constructor.

        Parameters:
        previous_results - Current results structure of the simulation run.
        job_server       - ParallelPython job server.
        """

        self.previous_results = previous_results
        self.job_server = job_server

        # will need the smec opd to mpd factor
        fts = self.previous_results['fts']
        self.smec_opd_to_mpd = fts['smec_opd_to_mpd']

        self.reduceint = []
        self.result = collections.OrderedDict()      

    def run(self):
#        print 'ReduceInterferogram.run'

        # namedtuple to hold UVspectrum result
        UVspectrum = collections.namedtuple(
          'uvspectrum', [
          'scan_number',
          'baseline_x',
          'baseline_y',
          'baseline_z',
          'interferogram',
          'opd', 
          'spectrum',
          'wavenumber'],
          verbose=False)

        # get observation list
        observe = self.previous_results['observe']
        obs_timeline = observe['observed_timeline']

        # Construct a dict holding lists of configurations for each
        # baseline observed
        baseline_configs = collections.defaultdict(list)
        observed_times = obs_timeline.keys()
        observed_times.sort()

        for t in observed_times:
            config = obs_timeline[t]
            b = (config.baseline_x, config.baseline_y, config.baseline_z)
            baseline_configs[b].append(config)

        baselines = baseline_configs.keys()
        baselines.sort()

        # for each baseline assemble config data into arrays
        baseline_scans = {}
        baseline_mean = {}
        baseline_uvspectrum = {}

        for b in baselines:
            data = []
            smec_position = []
            flag = []

            for config in baseline_configs[b]:
                data.append(config.data)
                smec_position.append(config.smec_position)
                flag.append(config.flag)

            data = np.array(data)
            smec_position = np.array(smec_position)
            flag = np.array(flag)        

            # identify scans within the data flow
            flag_sections = find_flagged_sections(flag)

            # extract scans from the data stream
            scans = []
            mark = 0
            for flag_section in flag_sections:
                interferogram = data[mark:flag_section[0]]
                opd = 100.0 * smec_position[mark:flag_section[0]] / \
                  self.smec_opd_to_mpd

                if len(interferogram):
                    axis = co.Axis(data=opd, title='OPD', units='cm')
                    scan = co.Spectrum(
                      data=interferogram,
                      flag=np.zeros(np.shape(interferogram), np.bool),
                      axis=axis,
                      title='scan interferogram',
                      units='')
                    scans.append(scan)

                mark = flag_section[1] + 1

            # deal with last scan
            interferogram = data[mark:]
            opd = 100.0 * smec_position[mark:] / self.smec_opd_to_mpd

            if len(interferogram):
                axis = co.Axis(data=opd, title='OPD', units='cm')
                scan = co.Spectrum(
                  data=interferogram,
                  flag=np.zeros(np.shape(interferogram), np.bool),
                  axis=axis, 
                  title='scan interferogram',
                  units='')
                scans.append(scan)
            
            baseline_scans[b] = scans

            # now obtain mean of scans for each baseline

            interferogram_sum = np.zeros(np.shape(scans[0].data))
            interferogram_n = np.zeros(np.shape(scans[0].data))
    
            for scan in scans:
                # assume opd ranges are the same but sort them into
                # ascending order
                opd = scan.axis.data
                opd_sort = np.argsort(opd)
                interferogram_sum += scan.data[opd_sort]
                interferogram_n += 1.0

            interferogram_mean = interferogram_sum / interferogram_n

            axis = co.Axis(data=opd[opd_sort], title='OPD', units='cm')
            scan_mean = co.Spectrum(
              data=interferogram_mean,
              flag=np.zeros(np.shape(interferogram_mean), np.bool),
              axis=axis, 
              title='mean interferogram',
              units='')

            baseline_mean[b] = scan_mean

            # now convert interferogram to spectrum at this u-v position
            # first, shift interferogram so that 0 opd is at index 0
            # must ask David N how to do this properly
            postzero = []
            prezero = []
            postzero_opd = []
            prezero_opd = []
            zero_opd_found = False

            for i, opd in enumerate(scan_mean.axis.data):
                if abs(opd) < 1e-10:
                    zero_opd_found = True
                if not zero_opd_found:
                    prezero.append(scan_mean.data[i])
                    prezero_opd.append(scan_mean.axis.data[i])
                else:
                    postzero.append(scan_mean.data[i])
                    postzero_opd.append(scan_mean.axis.data[i])
            
            shifted = postzero + prezero
            shifted_opd = postzero_opd + prezero_opd

            shifted = np.array(shifted)
            shifted_opd = np.array(shifted_opd)

            # spectrum is complex
            spectrum = np.fft.fft(shifted)
            wavenumber = np.fft.fftfreq(
              n=spectrum.size,
              d=abs(shifted_opd[1] - shifted_opd[0]))
      
            # save it
            uvspectrum = UVspectrum(
              None,
              b[0], 
              b[1],
              0.0,
              shifted,
              shifted_opd,
              spectrum,
              wavenumber)             

            baseline_uvspectrum[b] = uvspectrum

        self.result['baseline_scans'] = baseline_scans
        self.result['baseline_mean'] = baseline_mean
        self.result['baseline_uvspectrum'] = baseline_uvspectrum
        
        return self.result
                
    def __repr__(self):
        return '''
ReduceInterferogram : {num_uvspectra}
'''.format(
          num_uvspectra=len(self.result['baseline_uvspectrum']))

