import numpy as np

class Axis(object):
    """Class to hold axis information.
    """
    def __init__(self, data, title='', units=''):
        self.data = data.copy()
        self.title = title
        self.units = units


class Cube(object):
    """Class to hold 3d image, with axes, units, title, etc.
    """
    def __init__(self, data, axes=[], title='', units=''):
        self.data = data.copy()
        self.axes = axes
        self.title = title
        self.units = units


class Image(object):
    """Class to hold 2d image, with axes, units, title, etc.
    """
    def __init__(self, data, axes=[], title='', units='', info={}):
        self.data = data.copy()
        self.axes = axes
        self.title = title
        self.units = units
        self.info = info


class Spectrum(object):
    """Class to hold 1d spectrum, with axis, units, title, etc.
    """
    def __init__(self, data, flag=None, axis=None, title='', units=''):
        self.data = data.copy()
        if flag is None:
            self.flag = np.zeros(np.shape(data), np.bool)
        else:
            self.flag = flag.copy()
        self.axis = axis
        self.title = title
        self.units = units
