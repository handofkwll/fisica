from __future__ import absolute_import

import collections
import numpy as np
import matplotlib.pyplot as plt


import common.commonobjects as co

def find_flagged_sections(flag):
    """Routine to find ranges [start:end] of sections
    in the given array that are set true. Range means
    that points from flag[start] up to but not including
    flag[end] are True.

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
                flagsection[1] = i
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
        observe = self.previous_results['readfits']
        self.smec_opd_to_mpd = observe['smec_opd_to_mpd']
        obs_timeline = observe['observed_timeline']

        # Construct a dict holding lists of configurations for each
        # baseline observed
        baseline_configs = collections.defaultdict(list)
        observed_times = obs_timeline.keys()
        observed_times.sort()

        scanset = set()
        scans_data = {}
        data = None

        for t in observed_times:
            config = obs_timeline[t]
            scan = config.scan_number

            if scan not in scanset:
                if data is not None:
                    data = np.array(data)
                    smec_position = np.array(smec_position)
                    flag = np.array(flag)        
                    baseline_x = np.array(baseline_x)
                    baseline_y = np.array(baseline_y)
                    baseline_z = np.array(baseline_z)

                    smec_sorted = np.argsort(smec_position)
                    data = data[smec_sorted]
                    smec_position = smec_position[smec_sorted]
                    flag = flag[smec_sorted]
                    baseline_x = baseline_x[smec_sorted]
                    baseline_y = baseline_y[smec_sorted]
                    baseline_z = baseline_z[smec_sorted]

                    scans_data[current_scan] = (
                      data,
                      smec_position,
                      flag,
                      baseline_x,
                      baseline_y,
                      baseline_z)

                data = [config.data]
                smec_position = [config.smec_nominal_position]
                flag = [config.flag]
                baseline_x = [config.baseline_x]
                baseline_y = [config.baseline_y]
                baseline_z = [config.baseline_z]
                current_scan = scan
                scanset.update([scan])
            else:
                data.append(config.data)
                smec_position.append(config.smec_nominal_position)
                flag.append(config.flag)
                baseline_x.append(config.baseline_x)
                baseline_y.append(config.baseline_y)
                baseline_z.append(config.baseline_z)

        scans = scans_data.keys()
        scans.sort()
        scan_interferograms = {}
        scan_uvspectra = {}

        for scan in scans:
            data, smec_position, flag, baseline_x, baseline_y, \
              baseline_z = scans_data[scan]       

            interferogram = data[flag==False]
            opd = 100.0 * smec_position[flag==False] / self.smec_opd_to_mpd

            if len(interferogram):
                axis = co.Axis(data=opd, title='OPD', units='cm')
                scan_interferogram = co.Spectrum(
                  data=interferogram,
                  flag=np.zeros(np.shape(interferogram), np.bool),
                  axis=axis,
                  title='scan interferogram',
                  units='')
                scan_interferograms[scan] = scan_interferogram

                # now convert interferogram to spectrum at this u-v position
                # first, shift interferogram so that 0 opd is at index 0
                # must ask David N how to do this properly
                postzero = []
                prezero = []
                postzero_opd = []
                prezero_opd = []
                zero_opd_found = False

                for i, opd in enumerate(scan_interferogram.axis.data):
                    if abs(opd) < 1e-10:
                        zero_opd_found = True
                    if not zero_opd_found:
                        prezero.append(scan_interferogram.data[i])
                        prezero_opd.append(scan_interferogram.axis.data[i])
                    else:
                        postzero.append(scan_interferogram.data[i])
                        postzero_opd.append(scan_interferogram.axis.data[i])
            
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
                  scan_number=scan,
                  baseline_x = baseline_x,
                  baseline_y = baseline_y,
                  baseline_z = baseline_z,
                  interferogram=shifted,
                  opd=shifted_opd,
                  spectrum=spectrum,
                  wavenumber=wavenumber)             

                scan_uvspectra[scan] = uvspectrum

        self.result['scan_interferograms'] = scan_interferograms
        self.result['scan_uvspectra'] = scan_uvspectra
        
        return self.result
                
    def __repr__(self):
        return '''
ReduceInterferogram : {num_uvspectra}
'''.format(
          num_uvspectra=len(self.result['scan_uvspectra']))

