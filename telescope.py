"""This module contains the Telescope class.
"""

from __future__ import absolute_import

import collections
import numpy as np


class Telescope(object):
    """Class describing the flux collectors.

    Contains the methods:
    __init__
    run
    __repr__
    """

    def __init__(self, parameters):
        """Telescope constructor.

        Keyword arguments:
        parameters -- dict with parameters from the Excel configuration file.
        """
        self.parameters = parameters

        # read params from Telescope sheet
        telparm = parameters['substages']['Telescope']

        self.result = collections.OrderedDict()

        telparm = parameters['substages']['Telescope']
        self.result['m1_diameter'] = telparm['Primary mirror diameter']
        self.result['mirror_temp'] = telparm['Mirror temp']
        self.result['telescope_emissivity'] = telparm['Telescope emissivity (Etel)']
        self.result['n_telescope_mirrors'] = telparm['Number of mirrors in telescope design']
        self.result['c1_pointing_error_type'] = telparm['Collector 1 pointing error type']
        self.result['c2_pointing_error_type'] = 'Zero'
        self.result['c1_beam_model_type'] = telparm['Collector 1 beam model type']
        self.result['c2_beam_model_type'] = telparm['Collector 2 beam model type']
        self.result['beam_model_pol'] = telparm['Beam model polarization']
        self.result['field_rotator'] = bool(telparm['Field rotator'])

    def run(self):
        """Method returns a structure containing the derived Telescope
        parameters.
        """
        return self.result

    def __repr__(self):
        return '''
Telescope:
  m1 diam        : {diameter}
  mirror temp    : {temp}
  tel emissivity : {emissivity}
  n mirrors      : {nmirror}
  collector 1 pointing error : {pointing1}
  collector 2 pointing error : {pointing2}
  collector 1 beam model : {model1}
  collector 2 beam model : {model2}
  polarization           : {pol}
  field rotator          : {field_rot}
'''.format(
          diameter=self.result['m1_diameter'],
          temp=self.result['mirror_temp'],
          emissivity=self.result['telescope_emissivity'],
          nmirror=self.result['n_telescope_mirrors'],
          pointing1=self.result['c1_pointing_error_type'],
          pointing2=self.result['c2_pointing_error_type'],
          model1=self.result['c1_beam_model_type'],
          model2=self.result['c2_beam_model_type'],
          pol=self.result['beam_model_pol'],
          field_rot=self.result['field_rotator'])

