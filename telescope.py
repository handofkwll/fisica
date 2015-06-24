from __future__ import absolute_import

import collections
import numpy as np


class Telescope(object):
    """Class to describe the flux collectors.
    """

    def __init__(self, parameters):
        self.parameters = parameters

        # read params from Telescope sheet
        telparm = parameters['substages']['Telescope']

        self.result = collections.OrderedDict()

        telparm = parameters['substages']['Telescope']
        self.result['m1_diameter'] = telparm['Primary mirror diameter']
        self.result['mirror_temp'] = telparm['Mirror temp']
        self.result['telescope_emissivity'] = telparm['Telescope emissivity (Etel)']
        self.result['n_telescope_mirrors'] = telparm['Number of mirrors in telescope design']
        self.result['pointing_error_type'] = telparm['Pointing error type']
        self.result['beam_model_type'] = telparm['Beam model type']

    def run(self):
        #print 'Telescope.run'        
        return self.result

    def __repr__(self):
        return '''
Telescope:
  m1 diam        : {diameter}
  mirror temp    : {temp}
  tel emissivity : {emissivity}
  n mirrors      : {nmirror}
  pointing error : {pointing}
  beam model     : {model}
'''.format(
          diameter=self.result['m1_diameter'],
          temp=self.result['mirror_temp'],
          emissivity=self.result['telescope_emissivity'],
          nmirror=self.result['n_telescope_mirrors'],
          pointing=self.result['pointing_error_type'],
          model=self.result['beam_model_type'])

