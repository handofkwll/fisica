from __future__ import absolute_import

import collections
import numpy as np


class FTS(object):
    """Class to describe the Fourier Transform Spectrometer.
    """

    def __init__(self, parameters):
        self.parameters = parameters

        # read params from FTSpectrograph sheet
        ftsparm = parameters['substages']['FTSpectrograph']

        self.result = collections.OrderedDict()

        # slightly awkward getting values as they are keyed by 
        # the row number in the parent spreadsheet
        row = ftsparm['wn min [cm-1]'].keys()[0]
        self.result['wnmin'] = wnmin = ftsparm['wn min [cm-1]'][row]
        self.result['wnmax'] = wnmax = ftsparm['wn max [cm-1]'][row]
        self.result['specres'] = specres = ftsparm['Spectral Res'][row]

        # read params from FTSMechanical sheet
        ftsparm = parameters['substages']['FTSMechanical']
        row = ftsparm['V drive [cm/s]'].keys()[0]
        self.result['smec_velocity'] = smec_velocity = \
          ftsparm['V drive [cm/s]'][row]
        self.result['smec_verror_type'] = ftsparm['V error type'][row]
        self.result['nscans'] = int(ftsparm['Nscans'][row])
        self.result['smec_interscan_duration'] = \
          ftsparm['Interscan delay [s]'][row]

        # derived parameters
        self.result['delta_wn'] = delta_wn = wnmin / specres
        self.result['delta_opd'] = 1.0 / (2.0 * wnmax)
        # number of unaliased spectral points
        nspec = int(np.ceil(wnmax / delta_wn)) + 1
        # fts_nsample symmetric about 0 opd
        self.result['fts_nsample'] = fts_nsample = 2 * (nspec - 1)
        self.result['opd_max'] = opd_max = 1.0 / (2.0 * delta_wn)

        # unaliased spectral points sampled by FTS (wn=wavenumber[cm-1]
        i = np.arange(nspec)
        fts_wn = wnmax * i / float(nspec-1)

        # wn axis truncated at wnmin
        fts_wn_truncated = fts_wn[fts_wn >= wnmin]

        self.result['fts_wn'] = fts_wn
        self.result['fts_wn_truncated'] = fts_wn_truncated

        # smec parameters
        self.result['smec_opd_to_mpd'] = smec_opd_to_mpd = 4.0
        self.result['smec_start'] = smec_start = -opd_max / smec_opd_to_mpd
        self.result['smec_scan_duration'] = smec_scan_duration = \
          2.0 * abs(smec_start) / smec_velocity
    
    def run(self):
        #print 'FTS.run'        
        return self.result

    def __repr__(self):
        return '''
FTS:
  wnmin    : {wnmin}
  wnmax    : {wnmax}
  specres  : {specres}
  delta wn : {delta_wn}
  delta opd: {delta_opd}
  nsample  : {nsample}
  max opd  : {opd_max}
  v drive  : {vdrive}
  v error  : {verror}
  nscans   : {nscans}
  interscan: {interscan}
'''.format(
          wnmin=self.result['wnmin'],
          wnmax=self.result['wnmax'],
          specres=self.result['specres'],
          delta_wn=self.result['delta_wn'],
          delta_opd=self.result['delta_opd'],
          nsample=self.result['fts_nsample'],
          opd_max=self.result['opd_max'],
          vdrive=self.result['smec_velocity'],
          verror=self.result['smec_verror_type'],
          nscans=self.result['nscans'],
          interscan=self.result['smec_interscan_duration'])

